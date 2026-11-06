"""DeLong (generalised to two labels on the same patients, as in stats_delong_1b.py), paired bootstrap, Wilson CI."""
import numpy as np
from scipy.stats import norm
from sklearn.metrics import roc_auc_score


def influence(y, p):
    """AUC and per-patient influence terms z with sum(z_k * z_l) = DeLong covariance (S10/m + S01/n form)."""
    y, p = np.asarray(y).astype(int), np.asarray(p, dtype=float)
    pos, neg = p[y == 1], p[y == 0]
    m, n = len(pos), len(neg)
    cmp = (pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :])
    v10, v01 = cmp.mean(axis=1), cmp.mean(0)
    auc = v10.mean()
    z = np.empty(len(y))
    z[y == 1] = (v10 - auc) / np.sqrt(m * (m - 1))
    z[y == 0] = (v01 - auc) / np.sqrt(n * (n - 1))
    return auc, z


def delong(y_ref, p_ref, y_alt, p_alt):
    """AUC(alt) − AUC(ref) on the same patients (labels may differ): dict with both AUCs, delta, 95% CI, p"""
    a_ref, z_ref = influence(y_ref, p_ref)
    a_alt, z_alt = influence(y_alt, p_alt)
    se = np.sqrt(np.sum((z_alt - z_ref) ** 2))
    d = a_alt - a_ref
    return {"auc_ref": a_ref, "auc_alt": a_alt, "diff": d, "lo": d - 1.96 * se, "hi": d + 1.96 * se,
            "p": 2 * norm.sf(abs(d) / se) if se > 0 else 1.0}


def auc_ci(y, p):
    a, z = influence(y, p)
    se = np.sqrt(np.sum(z ** 2))
    return a, a - 1.96 * se, a + 1.96 * se


def paired_bootstrap(stat, n, b=2000, seed=0):
    """percentile CI of stat(idx) over patient resamples; stat returns a float (or nan to skip)"""
    r = np.random.default_rng(seed)
    v = np.array([stat(r.integers(0, n, n)) for _ in range(b)])
    v = v[np.isfinite(v)]
    return np.percentile(v, [2.5, 97.5])


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    ph = k / n
    c = (ph + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def safe_auc(y, p):
    return roc_auc_score(y, p) if len(np.unique(y)) == 2 else np.nan
