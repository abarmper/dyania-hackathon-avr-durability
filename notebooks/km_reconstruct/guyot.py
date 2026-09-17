"""Faithful Python port of Guyot et al. 2012 (BMC Med Res Methodol 12:9):
reconstructing individual patient data from a digitised Kaplan-Meier curve.

Inputs
  T, S        : digitised curve coordinates (T[0] = 0, S[0] = 1, S non-increasing, S in [0, 1])
  risk        : numbers at risk at times risk.times (must include 0); each time must be in T
  tot_events  : total number of events, if published
Information cases (chosen automatically):
  "all"            risk table with >= 2 times and tot_events
  "no_total"       risk table with >= 2 times, no tot_events
  "no_risk"        risk table at 0 only (n0) with tot_events
  "neither"        n0 only

Returns IPD as arrays (time, event) plus diagnostics.
"""
from __future__ import annotations

import numpy as np

from survival import km


class GuyotError(ValueError):
    pass


def _prepare(T, S, risk_times):
    T = np.asarray(T, float); S = np.asarray(S, float)
    if T[0] != 0:
        T = np.concatenate([[0.0], T]); S = np.concatenate([[1.0], S])
    # insert risk times not present, evaluating the step function (right-continuous)
    for rt in risk_times:
        if not np.any(np.isclose(T, rt)):
            idx = np.searchsorted(T, rt, side="right") - 1
            T = np.insert(T, idx + 1, rt); S = np.insert(S, idx + 1, S[idx])
    # enforce monotone non-increasing (pool adjacent violators, simple running min)
    S = np.minimum.accumulate(np.clip(S, 0.0, 1.0))
    return T, S


def reconstruct(T, S, risk_times, risk_n, tot_events=None, censor_at="coordinate", max_iter=300, obs_censor_times=None):
    """Run the Guyot algorithm. Returns (times, events, diag).

    obs_censor_times: digitised censor ("+") marks. When given, the last interval (which has no closing
    risk number) starts from the observed count and positions instead of Guyot's span-ratio guess;
    earlier intervals are unchanged because their risk numbers pin the censoring counts."""
    risk_times = list(risk_times); risk_n = [int(x) for x in risk_n]
    if len(risk_times) == 0 or risk_times[0] != 0:
        raise GuyotError("risk table must start at time 0 with n0")
    T, S = _prepare(T, S, risk_times)
    N = len(T)
    nint = len(risk_times)
    lower = [int(np.argmin(np.abs(T - rt))) for rt in risk_times]
    upper = [lower[i + 1] - 1 for i in range(nint - 1)] + [N - 1]

    if nint >= 2 and tot_events is not None:
        case = "all"
    elif nint >= 2:
        case = "no_total"
    elif tot_events is not None:
        case = "no_risk"
    else:
        case = "neither"

    d = np.zeros(N, int); cen = np.zeros(N, int); n = np.zeros(N + 1, int)
    n[0] = risk_n[0]
    ncensor = np.zeros(nint, int)
    placed = [[] for _ in range(nint)]
    warnings = []

    def run_interval(i, ncens, fixed=None):
        """Steps 2-3 for interval i with ncens censorings; fills d, cen, n in place. Returns n at upper_i+1.
        fixed: observed censor times to use first (sorted); any surplus is placed evenly."""
        lo, up = lower[i], upper[i]
        # step 2: even placement (or observed marks)
        t0, t1 = T[lo], (T[lower[i + 1]] if i + 1 < nint else T[up])
        cts = []
        if fixed:
            cts = list(sorted(fixed))[:ncens]
        extra = ncens - len(cts)
        if extra > 0 and t1 > t0:
            cts += [t0 + c * (t1 - t0) / (extra + 1) for c in range(1, extra + 1)]
        placed[i] = cts
        for k in range(lo, up + 1):
            cen[k] = 0
        for ct in cts:
            k = int(np.searchsorted(T, ct, side="right") - 1)
            k = min(max(k, lo), up)
            cen[k] += 1
        # step 3: events and numbers at risk
        n[lo] = risk_n[i]
        # running KM must continue from the previous interval's state
        S_run = 1.0
        for k in range(lo):
            if d[k] > 0 and n[k] > 0:
                S_run *= (1 - d[k] / n[k])
        for k in range(lo, up + 1):
            if n[k] <= 0:
                d[k] = 0; cen[k] = 0; n[k + 1] = 0; continue
            dk = int(round(n[k] * (1 - S[k] / S_run))) if S_run > 0 else 0
            dk = max(0, min(dk, n[k]))
            d[k] = dk
            if dk > 0:
                S_run *= (1 - dk / n[k])
            # censorings cannot exceed the patients left after the events at k
            cen[k] = int(min(cen[k], n[k] - dk))
            n[k + 1] = n[k] - dk - cen[k]
        return n[up + 1]

    # intervals with a following risk number (steps 1-5)
    for i in range(nint - 1):
        lo, lo_next = lower[i], lower[i + 1]
        ncens = int(round(S[lo_next] / S[lo] * risk_n[i])) - risk_n[i + 1] if S[lo] > 0 else 0
        ncens = max(0, ncens)
        best = None
        for it in range(max_iter):
            n_end = run_interval(i, ncens)
            gap = n_end - risk_n[i + 1]
            if best is None or abs(gap) < abs(best[0]):
                best = (gap, ncens)
            if gap == 0:
                break
            ncens = max(0, ncens + gap)
        else:
            warnings.append(f"interval {i}: risk-number match not exact (gap {best[0]}) after {max_iter} iterations")
            ncens = best[1]; run_interval(i, ncens)
        ncensor[i] = ncens

    # step 6: last interval
    i = nint - 1
    fixed_last = None
    if obs_censor_times is not None:
        fixed_last = [float(c) for c in obs_censor_times if c >= T[lower[i]] - 1e-9]
        ncens = min(len(fixed_last), risk_n[i])
        last_interval_rule = f"observed marks ({len(fixed_last)})"
    elif case in ("no_risk", "neither"):
        ncens = 0; last_interval_rule = "none (no risk table)"
    else:
        span_last = T[upper[i]] - T[lower[i]]
        span_prev = T[upper[i - 1]] - T[lower[0]] if i >= 1 else 1.0
        ncens = int(round(span_last / span_prev * ncensor[:i].sum())) if span_prev > 0 else 0
        ncens = min(max(ncens, 0), risk_n[i]); last_interval_rule = "Guyot span ratio"
    run_interval(i, ncens, fixed_last)
    ncensor[i] = ncens

    # steps 7-8: total events
    if tot_events is not None:
        events_before = int(d[:lower[i]].sum())
        tail_drops = bool(S[upper[i]] < S[lower[i]] - 1e-9)
        if events_before >= tot_events and case == "all":
            if tail_drops:
                # The published total is inconsistent with the digitised curve (which still falls in the
                # last interval). Zeroing the tail, as step 7 prescribes, would discard real curve shape;
                # keep the curve and report the inconsistency instead.
                warnings.append(f"events before last interval ({events_before}) >= published total ({tot_events}) "
                                f"although the curve still falls after t={T[lower[i]]:.2f}; total-events constraint ignored")
            else:
                for k in range(lower[i], upper[i] + 1):
                    d[k] = 0; cen[k] = 0; n[k + 1] = n[lower[i]]
                if events_before > tot_events:
                    warnings.append(f"events before last interval ({events_before}) exceed total ({tot_events})")
        else:
            for it in range(max_iter):
                total = int(d.sum())
                if total == tot_events:
                    break
                ncens_new = max(0, ncens + (total - tot_events))
                if ncens_new == ncens and total < tot_events:
                    warnings.append(f"cannot reach total events {tot_events}; got {total} with no censoring left")
                    break
                if ncens_new == ncens:
                    # try the other direction once
                    ncens_new = max(0, ncens + (total - tot_events))
                ncens = ncens_new
                run_interval(i, ncens, fixed_last)
            ncensor[i] = ncens

    # build IPD
    times, events = [], []
    for k in range(N):
        times += [T[k]] * int(d[k]); events += [1] * int(d[k])
        if censor_at == "placed":
            pass
        else:
            times += [T[k]] * int(cen[k]); events += [0] * int(cen[k])
    if censor_at == "placed":
        for i in range(nint):
            for ct in placed[i]:
                times.append(ct); events.append(0)
        # placed censorings may be fewer than cen if rounding; top up at coordinates
        extra = int(cen.sum()) - sum(len(p) for p in placed)
        if extra > 0:
            for k in range(N):
                while extra > 0 and cen[k] > 0:
                    times.append(T[k]); events.append(0); extra -= 1; cen[k] -= 1
    times = np.array(times, float); events = np.array(events, int)
    # patients still at risk at the end are censored at the last coordinate
    remaining = risk_n[0] - len(times)
    if remaining > 0:
        times = np.concatenate([times, np.full(remaining, T[-1])]); events = np.concatenate([events, np.zeros(remaining, int)])
    elif remaining < 0:
        warnings.append(f"IPD count {len(times)} exceeds n0 {risk_n[0]}")

    diag = {
        "case": case, "n_ipd": int(len(times)), "total_events": int(events.sum()),
        "ncensor_per_interval": ncensor.tolist(), "warnings": warnings,
        "risk_times": risk_times, "risk_n": risk_n,
        "n_at_risk_reconstructed": [int(n[lower[i]]) for i in range(nint)],
        "last_interval_rule": last_interval_rule,
        "T": T.tolist(), "S": S.tolist(), "d": d.tolist(), "cen": cen.tolist(),
    }
    return times, events, diag


def reconstructed_curve(times, events):
    return km(times, events)


def exact_inversion(T, S, n0, t_end=None, max_multiplicity=6):
    """Vector-only reconstruction: every KM step is an exact ratio S_k/S_{k-1} = 1 - d_k/n_k.

    For each drop choose the smallest event multiplicity m such that the implied risk set
    n = m / (1 - S_k/S_{k-1}) does not exceed the risk set left after the previous step;
    the difference between consecutive risk sets is censoring, placed evenly between the steps.
    Requires vector-precision coordinates (line-width errors of ~0.7 pp break it late in follow-up,
    so the tail with n < 20 is approximate). Returns (times, events, diag)."""
    T = np.asarray(T, float); S = np.asarray(S, float)
    times, events = [], []
    n_prev = int(n0); S_prev = 1.0; t_prev = 0.0
    n_trace = [(0.0, n_prev)]; warnings = []
    for k in range(1, len(T)):
        if S[k] >= S_prev - 1e-9:
            continue
        rel = 1.0 - S[k] / S_prev
        m = 1
        for mm in range(1, max_multiplicity + 1):
            if mm / rel <= n_prev + 0.5:
                m = mm; break
        else:
            m = 1; warnings.append(f"t={T[k]:.2f}: drop {rel*100:.2f}% too small for n<={n_prev}; using 1 event")
        n_k = int(min(round(m / rel), n_prev))
        cens = n_prev - n_k
        if cens > 0:
            for c in range(1, cens + 1):
                times.append(t_prev + c * (T[k] - t_prev) / (cens + 1)); events.append(0)
        times += [T[k]] * m; events += [1] * m
        n_prev = n_k - m; S_prev = S[k]; t_prev = T[k]
        n_trace.append((float(T[k]), n_k))
    end = float(t_end if t_end is not None else T[-1])
    times += [end] * n_prev; events += [0] * n_prev
    diag = {"case": "exact_inversion", "n_ipd": len(times), "total_events": int(sum(events)),
            "n_at_risk_trace": n_trace, "warnings": warnings}
    return np.array(times, float), np.array(events, int), diag
