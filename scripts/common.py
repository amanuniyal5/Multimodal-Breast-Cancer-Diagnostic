"""Shared helpers: load TITAN slide embeddings, the linear-probe classifier, bootstrap CIs."""
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def load_titan(path):
    """Patient-level embedding = mean of that patient's DX-slide TITAN embeddings (768-d)."""
    t = pq.read_table(path, columns=["case_submitter_id", "svs_filename", "embeddings"]).to_pandas()
    t = t[t.svs_filename.str.contains("-DX")]
    X = np.stack([np.asarray(list(e[0]), dtype=np.float32) for e in t.embeddings])
    return pd.DataFrame(X, index=t.case_submitter_id).groupby(level=0).mean()


def linear_probe():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(C=0.01, max_iter=3000, class_weight="balanced"))


def bootstrap_auc_ci(y, p, n=1000, seed=0):
    r = np.random.default_rng(seed)
    aucs = []
    for _ in range(n):
        i = r.integers(0, len(y), len(y))
        if len(set(y[i])) < 2:
            continue
        aucs.append(roc_auc_score(y[i], p[i]))
    return np.percentile(aucs, [2.5, 97.5])
