"""Round trip: simulate IPD -> KM coordinates + risk table + total events -> Guyot -> compare."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from guyot import reconstruct  # noqa: E402
from survival import km, step_eval  # noqa: E402


def _simulate(seed=0, n=300):
    rng = np.random.default_rng(seed)
    T = rng.weibull(1.7, n) * 9.0
    C = rng.uniform(2.0, 12.0, n)
    t = np.minimum(T, C); e = (T <= C).astype(int)
    return t, e


def _risk_table(t, times):
    return [int((t >= rt).sum()) for rt in times]


@pytest.mark.parametrize("case", ["all", "no_total", "no_risk", "neither"])
def test_roundtrip(case):
    t, e = _simulate(seed=3)
    kt, kS = km(t, e)
    risk_times = [0, 2, 4, 6, 8, 10] if case in ("all", "no_total") else [0]
    risk_n = _risk_table(t, risk_times)
    tot = int(e.sum()) if case in ("all", "no_risk") else None
    rt, re_, diag = reconstruct(kt, kS, risk_times, risk_n, tot_events=tot)
    assert diag["case"] == case
    rkt, rkS = km(rt, re_)
    grid = np.linspace(0, 10, 41)
    diff = np.max(np.abs(step_eval(grid, kt, kS) - step_eval(grid, rkt, rkS)))
    tol = {"all": 0.02, "no_total": 0.03, "no_risk": 0.05, "neither": 0.08}[case]
    assert diff <= tol, f"max KM diff {diff:.4f} > {tol} ({case})"
    if tot is not None:
        assert diag["total_events"] == tot, diag
    if case in ("all", "no_total"):
        assert diag["n_at_risk_reconstructed"] == risk_n, diag["n_at_risk_reconstructed"]
    assert diag["n_ipd"] == len(t)


def test_survival_probabilities_close_all_information():
    t, e = _simulate(seed=11, n=500)
    kt, kS = km(t, e)
    risk_times = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    rt, re_, diag = reconstruct(kt, kS, risk_times, _risk_table(t, risk_times), tot_events=int(e.sum()))
    rkt, rkS = km(rt, re_)
    for tt in (1, 3, 5, 7, 9):
        a = step_eval(tt, kt, kS)[0]; b = step_eval(tt, rkt, rkS)[0]
        assert abs(a - b) <= 0.01, (tt, a, b)
