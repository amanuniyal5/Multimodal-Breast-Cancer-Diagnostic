"""Outer folds (fixed, data/splits/folds.csv), inner folds derived from the same file, label-efficiency subsets.

Inner folds: for an outer split at seed s, the training patients are re-split by their fold_s{(s+1) mod 5} column, so
inner folds are site-grouped, deterministic and library-independent (no new random splits).
"""
import numpy as np
import pandas as pd

from common import site_folds

_F = {}


def _table(path="data/splits/folds.csv"):
    if path not in _F:
        _F[path] = pd.read_csv(path, index_col="patient")
    return _F[path]


def outer(patients, seed):
    return site_folds(list(patients), seed)


def inner(train_patients, seed):
    """[(tr_idx, va_idx)] positional within train_patients; groups with no training patients are dropped"""
    f = _table().loc[list(train_patients), f"fold_s{(seed + 1) % 5}"].to_numpy()
    return [(np.where(f != k)[0], np.where(f == k)[0]) for k in sorted(set(f)) if 0 < (f == k).sum() < len(f)]


def subsample(y_train, n, rep, seed):
    """stratified random subset of n training positions (all if n >= len), reproducible per (seed, rep)"""
    y_train = np.asarray(y_train)
    if n >= len(y_train):
        return np.arange(len(y_train))
    r = np.random.default_rng(1000 * seed + rep)
    out = []
    for c in np.unique(y_train):
        idx = np.where(y_train == c)[0]
        k = max(2, int(round(n * len(idx) / len(y_train))))
        out.append(r.choice(idx, min(k, len(idx)), replace=False))
    return np.sort(np.concatenate(out))
