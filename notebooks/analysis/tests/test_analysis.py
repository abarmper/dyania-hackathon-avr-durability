import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import analysis_config as C  # noqa: E402
from data import baseline_frame, echo_wide, history_arrays, load_tables  # noqa: E402
from features import build_landmark_table, history_features  # noqa: E402
from models import DiscreteTimeHazard  # noqa: E402
from outcomes import known_times, make_instances, person_periods  # noqa: E402


@pytest.fixture(scope="module")
def cohort():
    T = load_tables(); w = echo_wide(T["echo_visits"]); h = history_arrays(w, T["patients_valves"].valve_id.to_numpy())
    L = build_landmark_table(h, T["patients_valves"])
    return T, h, L


def test_online_equals_batch(cohort):
    T, h, L = cohort
    base = baseline_frame(T["patients_valves"]).set_index("valve_id").loc[h.valve_ids].reset_index()
    rows = np.where(h.n_echoes >= 3)[0][:200]
    online = history_features(h, 3, base, rows)
    batch = L[(L.k == 3) & L.valve_id.isin(online.valve_id)].sort_values("valve_id").reset_index(drop=True)
    online = online.sort_values("valve_id").reset_index(drop=True)
    for c in ("ewma_log_mpg", "slope_since_ref", "prior_positive_count", "consecutive_worsening", "delta_mpg"):
        np.testing.assert_allclose(online[c].to_numpy(float), batch[c].to_numpy(float), equal_nan=True)


def test_outcome_timing(cohort):
    T, h, L = cohort
    kt = known_times(T["labels"], events=T["events"])
    assert (kt.cause == 1).sum() >= 1641          # confirmed labels plus SVD-driven reinterventions that truncated follow-up
    inst = make_instances(L, T["labels"], events=T["events"])
    assert (inst.time_to > 0).all()
    assert (inst.T_known > inst.landmark_t).all()
    X, y, rep = person_periods(inst.head(2000), C.FULL)
    assert set(np.unique(y)) <= {0, 1, 2}
    per = X.groupby(rep)["period"].max().to_numpy(); tt = inst.head(2000).time_to.to_numpy()
    assert (per == np.minimum(np.ceil(tt / C.PERIOD_Y - 1e-9), 10)).all()


def test_cif_monotone_bounded(cohort):
    T, h, L = cohort
    inst = make_instances(L, T["labels"], events=T["events"]).sample(4000, random_state=0).reset_index(drop=True)
    m = DiscreteTimeHazard(C.BASELINE + ["landmark_t", "mpg_now", "delta_mpg"], max_iter=60).fit(inst)
    cif = m.predict_cif(inst.head(500), [1, 2, 3, 4, 5])
    assert (np.diff(cif, axis=1) >= -1e-12).all() and (cif <= 1).all() and (cif >= 0).all()
    cif2 = m.predict_cif(inst.head(500), [5], cause=2)
    assert ((cif[:, -1] + cif2[:, 0]) <= 1 + 1e-9).all()
