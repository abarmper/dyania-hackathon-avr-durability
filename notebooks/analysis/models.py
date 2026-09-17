"""Models. DiscreteTimeHazard = boosted multinomial cause-specific hazard on person-period rows (C2a);
CoxCauseSpecific = cause-specific Cox + cumulative incidence on the implant-time instances (C2b, HR table);
TimeSinceImplantOnly and RiskFactorCoxB4 = the two comparators of docs/03 section 0."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier

import analysis_config as C
from outcomes import person_periods


class DiscreteTimeHazard:
    """P(cause c in period j | alive at start of j, features at the landmark). Softmax over {none, SVD, competing}, so
    h1 + h2 <= 1 and the CIF (sum_j S_{j-1} h1_j) is monotone and bounded by construction."""

    def __init__(self, feature_cols, horizon_y=C.LANDMARK_HORIZON_Y, period_y=C.PERIOD_Y, max_iter=600, learning_rate=0.05,
                 max_leaf_nodes=31, min_samples_leaf=100, l2=1.0, seed=C.SEED):
        self.feature_cols = list(feature_cols); self.horizon_y = horizon_y; self.period_y = period_y
        self.J = int(round(horizon_y / period_y))
        self.clf = HistGradientBoostingClassifier(loss="log_loss", max_iter=max_iter, learning_rate=learning_rate, max_leaf_nodes=max_leaf_nodes,
                                                  min_samples_leaf=min_samples_leaf, l2_regularization=l2, early_stopping=False,
                                                  categorical_features=None, random_state=seed)

    @property
    def cols(self):
        return self.feature_cols + ["period"]

    def fit(self, inst_train: pd.DataFrame, inst_val: pd.DataFrame | None = None):
        X, y, _ = person_periods(inst_train, self.feature_cols, self.horizon_y, self.period_y)
        if inst_val is not None and len(inst_val):
            Xv, yv, _ = person_periods(inst_val, self.feature_cols, self.horizon_y, self.period_y)
            self.clf.set_params(early_stopping=True, validation_fraction=None, n_iter_no_change=30, tol=1e-5)
            self.clf.fit(X[self.cols].to_numpy(float), y, X_val=Xv[self.cols].to_numpy(float), y_val=yv)
        else:
            self.clf.fit(X[self.cols].to_numpy(float), y)
        self.n_train_rows_ = len(X); self.n_iter_ = getattr(self.clf, "n_iter_", None)
        return self

    def predict_hazards(self, inst: pd.DataFrame, n_periods: int | None = None) -> np.ndarray:
        """(n, J, 3) probabilities [none, SVD, competing] per period."""
        J = self.J if n_periods is None else min(n_periods, self.J)
        n = len(inst)
        Xb = inst[self.feature_cols].to_numpy(float)
        X = np.repeat(Xb, J, axis=0)
        period = np.tile(np.arange(1, J + 1), n)[:, None].astype(float)
        P = self.clf.predict_proba(np.hstack([X, period]))
        out = np.zeros((n * J, 3)); out[:, self.clf.classes_.astype(int)] = P
        return out.reshape(n, J, 3)

    def predict_cif(self, inst: pd.DataFrame, deltas_y, cause: int = 1) -> np.ndarray:
        """Cumulative incidence of `cause` by each delta (years) after the landmark; (n, len(deltas))."""
        H = self.predict_hazards(inst)
        h1, h2 = H[:, :, 1], H[:, :, 2]
        S_prev = np.cumprod(np.hstack([np.ones((len(inst), 1)), 1.0 - h1 - h2]), axis=1)[:, :-1]
        inc = S_prev * (h1 if cause == 1 else h2)
        cif = np.cumsum(inc, axis=1)
        idx = np.clip(np.round(np.asarray(deltas_y, float) / self.period_y).astype(int), 1, self.J) - 1
        return cif[:, idx]


class TimeSinceImplantOnly(DiscreteTimeHazard):
    """Null comparator: hazard depends on time since implant (landmark time and period) only."""

    def __init__(self, **kw):
        super().__init__(feature_cols=["landmark_t"], **kw)


class CoxCauseSpecific:
    """Cause-specific Cox models (lifelines) on implant-time instances; CIF_1(t) = int S(u-) dH_1(u) with
    S = exp(-(H_1 + H_2)); gives the hazard-ratio table (D5) and the published-risk-factor comparator (B4 set)."""

    def __init__(self, feature_cols, penalizer=0.01):
        self.feature_cols = list(feature_cols); self.penalizer = penalizer

    def _design(self, inst):
        X = inst[self.feature_cols].copy()
        X = X.fillna(X.median(numeric_only=True))
        return X

    def fit(self, inst: pd.DataFrame):
        from lifelines import CoxPHFitter
        X = self._design(inst)
        self.mu_ = X.mean(); self.sd_ = X.std().replace(0, 1.0)
        Z = (X - self.mu_) / self.sd_
        self.models_ = {}
        for cause in (1, 2):
            d = Z.copy(); d["T"] = inst.time_to.to_numpy(); d["E"] = (inst.cause.to_numpy() == cause).astype(int)
            m = CoxPHFitter(penalizer=self.penalizer).fit(d, "T", "E")
            self.models_[cause] = m
        return self

    def hr_table(self) -> pd.DataFrame:
        m = self.models_[1].summary
        out = pd.DataFrame({"feature": m.index, "HR_per_SD": np.exp(m["coef"]), "CI_low": np.exp(m["coef lower 95%"]), "CI_high": np.exp(m["coef upper 95%"]), "p": m["p"]})
        return out.sort_values("p").reset_index(drop=True)

    def predict_cif(self, inst: pd.DataFrame, taus, cause: int = 1) -> np.ndarray:
        Z = (self._design(inst) - self.mu_) / self.sd_
        H0 = {c: self.models_[c].baseline_cumulative_hazard_ for c in (1, 2)}
        times = np.union1d(H0[1].index.to_numpy(), H0[2].index.to_numpy())
        H01 = np.interp(times, H0[1].index.to_numpy(), H0[1].iloc[:, 0].to_numpy(), left=0.0)
        H02 = np.interp(times, H0[2].index.to_numpy(), H0[2].iloc[:, 0].to_numpy(), left=0.0)
        lp1 = self.models_[1].predict_partial_hazard(Z).to_numpy(); lp2 = self.models_[2].predict_partial_hazard(Z).to_numpy()
        dH01 = np.diff(np.r_[0.0, H01]); dH02 = np.diff(np.r_[0.0, H02])
        out = np.zeros((len(inst), len(taus)))
        for i_t, tau in enumerate(taus):
            m = times <= tau + 1e-9
            H1 = np.outer(lp1, np.cumsum(dH01[m])); H2 = np.outer(lp2, np.cumsum(dH02[m]))
            S_prev = np.exp(-(H1 + H2)); S_prev = np.hstack([np.ones((len(inst), 1)), S_prev[:, :-1]])
            inc = S_prev * np.outer(lp1 if cause == 1 else lp2, dH01[m] if cause == 1 else dH02[m])
            out[:, i_t] = inc.sum(1)
        return np.clip(out, 0, 1)


class RiskFactorCoxB4(CoxCauseSpecific):
    def __init__(self, penalizer=0.01):
        super().__init__(C.B4_SET, penalizer)
