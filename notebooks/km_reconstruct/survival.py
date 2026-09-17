"""Small survival toolkit: Kaplan-Meier, Aalen-Johansen, step evaluation, parametric MLE.

Pure numpy/scipy so it runs anywhere; lifelines is used only for cross-checks in tests.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def km(times, events):
    """Kaplan-Meier estimate. Returns (t, S) with t[0]=0, S[0]=1, one row per distinct event time."""
    t = np.asarray(times, float); e = np.asarray(events, int)
    order = np.argsort(t, kind="stable"); t = t[order]; e = e[order]
    ut = np.unique(t[e == 1])
    S = 1.0; ts = [0.0]; Ss = [1.0]
    n_at_risk = len(t)
    idx = 0
    for u in ut:
        # number at risk just before u: subjects with time >= u
        n = int((t >= u).sum())
        d = int(((t == u) & (e == 1)).sum())
        if n > 0:
            S *= (1 - d / n)
        ts.append(float(u)); Ss.append(float(S))
    return np.array(ts), np.array(Ss)


def aj(times, states):
    """Aalen-Johansen cumulative incidence. states: 0 censored, 1..K event types.
    Returns dict {k: (t, CIF_k)} with t[0]=0, CIF[0]=0."""
    t = np.asarray(times, float); s = np.asarray(states, int)
    order = np.argsort(t, kind="stable"); t = t[order]; s = s[order]
    ut = np.unique(t[s > 0])
    ks = sorted(int(k) for k in np.unique(s[s > 0]))
    S_prev = 1.0
    out = {k: ([0.0], [0.0]) for k in ks}
    for u in ut:
        n = int((t >= u).sum())
        if n == 0:
            continue
        d_all = int(((t == u) & (s > 0)).sum())
        for k in ks:
            d_k = int(((t == u) & (s == k)).sum())
            out[k][0].append(float(u))
            out[k][1].append(out[k][1][-1] + S_prev * d_k / n)
        S_prev *= (1 - d_all / n)
    return {k: (np.array(v[0]), np.array(v[1])) for k, v in out.items()}


def step_eval(t_grid, t, y, side="right"):
    """Evaluate a right-continuous step function (t, y) on t_grid.
    side='right': value at t includes jumps at t; side='left': value just before t."""
    t = np.asarray(t, float); y = np.asarray(y, float); tg = np.atleast_1d(np.asarray(t_grid, float))
    idx = np.searchsorted(t, tg, side="right" if side == "right" else "left") - 1
    idx = np.clip(idx, 0, len(y) - 1)
    return y[idx]


def _nll_weibull(params, t, e):
    log_lam, log_k = params
    lam, k = np.exp(log_lam), np.exp(log_k)
    z = t / lam
    ll = np.sum(e * (np.log(k / lam) + (k - 1) * np.log(np.maximum(z, 1e-300)))) - np.sum(z ** k)
    return -ll


def fit_weibull(times, events):
    """MLE for Weibull S(t) = exp(-(t/lambda)^k). Returns dict(scale=lambda, shape=k, nll, aic)."""
    t = np.asarray(times, float); e = np.asarray(events, int)
    t = np.maximum(t, 1e-6)
    x0 = np.array([np.log(max(np.mean(t), 1e-3) * 1.5), 0.0])
    res = minimize(_nll_weibull, x0, args=(t, e), method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-8, "fatol": 1e-10})
    lam, k = np.exp(res.x)
    return {"scale": float(lam), "shape": float(k), "nll": float(res.fun), "aic": float(2 * res.fun + 4)}


def _nll_loglogistic(params, t, e):
    log_a, log_b = params
    a, b = np.exp(log_a), np.exp(log_b)
    z = (t / a) ** b
    S = 1 / (1 + z)
    # hazard h = (b/t) z / (1+z)
    log_h = np.log(b) - np.log(t) + np.log(np.maximum(z, 1e-300)) - np.log1p(z)
    ll = np.sum(e * log_h) + np.sum(np.log(S))
    return -ll


def fit_loglogistic(times, events):
    """MLE for log-logistic S(t) = 1 / (1 + (t/alpha)^beta)."""
    t = np.asarray(times, float); e = np.asarray(events, int)
    t = np.maximum(t, 1e-6)
    x0 = np.array([np.log(max(np.median(t), 1e-3) * 1.5), 0.5])
    res = minimize(_nll_loglogistic, x0, args=(t, e), method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-8, "fatol": 1e-10})
    a, b = np.exp(res.x)
    return {"alpha": float(a), "beta": float(b), "nll": float(res.fun), "aic": float(2 * res.fun + 4)}


def weibull_S(t, scale, shape):
    return np.exp(-(np.asarray(t, float) / scale) ** shape)


def loglogistic_S(t, alpha, beta):
    return 1.0 / (1.0 + (np.asarray(t, float) / alpha) ** beta)
