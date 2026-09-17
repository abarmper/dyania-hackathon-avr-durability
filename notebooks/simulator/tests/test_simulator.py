import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "km_reconstruct"))

from labels import capodanno, confirmed, dvir, first_time, varc3_stage  # noqa: E402
from params import default_params  # noqa: E402
from physics import eoa_from_mpg, fit_physics_constants, mpg_from_eoa, ppm_class  # noqa: E402
from simulate import simulate_cohort  # noqa: E402


@pytest.fixture(scope="module")
def params():
    p = default_params(); fit_physics_constants(p); return p


def test_physics_round_trip(params):
    eoa = np.array([0.8, 1.2, 1.6, 2.2]); Q = np.array([250.0, 280.0, 300.0, 350.0]); c = 3.2e-4
    m = mpg_from_eoa(eoa, Q, c); back = eoa_from_mpg(m, Q, c)
    assert np.allclose(back, eoa, rtol=1e-9)
    assert np.all(np.diff(m) < 0)  # bigger orifice, lower gradient


def test_ppm_boundaries():
    assert list(ppm_class([0.86, 0.85, 0.66, 0.65, 0.64], [25] * 5)) == [0, 1, 1, 2, 2]
    assert list(ppm_class([0.71, 0.70, 0.56, 0.55], [32] * 4)) == [0, 1, 1, 2]


def test_ph_weibull_sampler():
    rng = np.random.default_rng(0); k, lam, n = 2.5, 10.0, 200_000
    for hr in (0.5, 1.0, 2.0):
        T = lam * (-np.log(rng.random(n)) / hr) ** (1 / k)
        t = np.array([3.0, 6.0, 9.0]); emp = np.array([(T > x).mean() for x in t]); theo = np.exp(-hr * (t / lam) ** k)
        assert np.max(np.abs(emp - theo)) < 0.004


def test_label_rules_hand_built():
    # rows: flicker (labelled then absent), confirmed, missing channels, regurgitant
    mpg = np.array([[24.0, 12.0, 13.0], [25.0, 27.0, np.nan], [24.0, np.nan, np.nan], [12.0, 12.0, 13.0]])
    eoa = np.array([[1.0, 1.5, 1.5], [1.0, 0.9, np.nan], [np.nan, np.nan, np.nan], [1.5, 1.5, 1.5]])
    dvi = np.array([[0.30, 0.45, 0.45], [0.30, 0.28, np.nan], [np.nan, np.nan, np.nan], [0.45, 0.45, 0.45]])
    ar = np.array([[0, 0, 0], [0, 0, np.nan], [0, np.nan, np.nan], [2, 2, 3]], float)
    ref = dict(mpg=np.full((4, 1), 12.0), eoa=np.full((4, 1), 1.5), dvi=np.full((4, 1), 0.45), ar=np.zeros((4, 1)))
    st, sth = varc3_stage(mpg, eoa, dvi, ar, ref["mpg"], ref["eoa"], ref["dvi"], ref["ar"])
    assert st[0, 0] == 2 and st[0, 1] == 0                       # flicker
    assert st[1, 0] == 2 and st[1, 1] == 2                       # persists
    assert np.isnan(st[2, 0]) and sth[2, 0] == 2                 # full definition NA without EOA/DVI, haemodynamic variant fires
    assert st[3, 0] == 2 and st[3, 2] == 3                       # regurgitant mode
    t = np.array([[1, 2, 3]] * 4, float)
    c2 = confirmed(st, t, 2)
    assert not c2[0].any() and c2[1, 0] and not c2[2].any()
    assert first_time(c2, t)[1] == 1.0 and np.isinf(first_time(c2, t)[0])
    assert capodanno(mpg, ref["mpg"])[0, 0] == 2 and dvir(mpg, eoa, dvi, ref["mpg"], ref["eoa"], ref["dvi"])[1, 1] == 2


def test_cohort_consistency_and_reproducibility(params):
    t0 = time.time(); a = simulate_cohort(params, n=3000, seed=7); dt = time.time() - t0
    b = simulate_cohort(params, n=3000, seed=7)
    assert dt < 30
    pd.testing.assert_frame_equal(a["patients_valves"], b["patients_valves"])
    pd.testing.assert_frame_equal(a["echo_visits"], b["echo_visits"])
    pv, ev, events, lab = a["patients_valves"], a["echo_visits"], a["events"], a["labels"]
    end = pv.set_index("valve_id").followup_end_years
    assert (ev.visit_time_years <= end.loc[ev.valve_id].to_numpy() + 1e-9).all()      # no echo after follow-up end
    dead = events[events.event_type == "death"].set_index("valve_id").time_years
    assert np.allclose(end.loc[dead.index], dead, atol=1e-6)                             # death ends follow-up
    assert (lab.primary_time <= lab.followup_end_years + 1e-9).all()
    assert set(ev.measurement) >= {"av_mean_gradient_mmHg", "av_area_cm2", "doppler_velocity_index"}
    assert pv.ref_mean_gradient_mmHg.between(1, 70).all()
    assert pv.ref_mean_gradient_mmHg.between(4, 40).mean() > 0.95
