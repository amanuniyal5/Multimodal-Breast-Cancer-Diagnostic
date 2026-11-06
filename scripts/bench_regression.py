"""X3: regression on mRNA expression from H&E (docs/DESIGN_SPEC.md §5.3), mean-pool features + ridge.

Targets log2 ESR1 / PGR / ERBB2. Ridge alpha chosen on inner folds (bench/folds.inner), 5 seeds, predictions averaged.
  regression.csv   per encoder x gene: Pearson r, Spearman rho, R² (out-of-fold) with patient-bootstrap 95% CIs
  reg_to_clf.csv   predicted expression used as a score: AUC vs the IHC and molecular labels, compared (DeLong, same
                   patients) with the classifier trained on that label (X1 oof.csv.gz); measured mRNA as the ceiling
  oof_reg.csv.gz   seed-averaged predictions
Usage: python scripts/bench_regression.py --x1 results/bench/x1 --out results/bench/x3 [--encoders ...]
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import folds, registry, stats, tasks  # noqa: E402
from bench.metrics import regression  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--encoders", default="")
ap.add_argument("--x1", default="results/bench/x1")
ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--alphas", default="10,100,1000,10000,100000")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--n-jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
ap.add_argument("--out", default="results/bench/x3")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()
ALPHAS = [float(x) for x in a.alphas.split(",")]
GENES = {"ESR1_log2": [("ER_IHC", "ESR1_RNApos")], "PGR_log2": [("PR_IHC", "PGR_prevmatch")],
         "ERBB2_log2": [("HER2_STRICT", "ERBB2_AMP")]}


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


L = tasks.labels()
ENCS = a.encoders.split(",") if a.encoders else [e for e in registry.available() if e != "resnet50"]
FE = {e: registry.load_means(e) for e in ENCS}
common = set(L.index) & set(pd.read_csv("data/splits/folds.csv").patient)
for e in ENCS:
    common &= set(FE[e].index)
common = sorted(common)
log(f"encoders {ENCS}; cohort {len(common)}")


def ridge(alpha):
    return make_pipeline(StandardScaler(), Ridge(alpha=alpha))


def fit_predict(X, y, pts, seed):
    p = np.zeros(len(y))
    for tr, te in folds.outer(pts, seed):
        inner = folds.inner(pts[tr], seed)
        err = []
        for al in ALPHAS:  # inner-CV mean squared error per alpha
            e = 0.0
            for itr, iva in inner:
                m = ridge(al).fit(X[tr][itr], y[tr][itr])
                e += ((m.predict(X[tr][iva]) - y[tr][iva]) ** 2).sum()
            err.append(e)
        p[te] = ridge(ALPHAS[int(np.argmin(err))]).fit(X[tr], y[tr]).predict(X[te])
    return p


def run(enc, gene):
    d = L.loc[common, [gene]].dropna()
    pts, y = d.index.to_numpy(), d[gene].to_numpy(dtype=float)
    X = FE[enc].loc[pts].values
    return enc, gene, pts, y, np.mean([fit_predict(X, y, pts, s) for s in range(a.seeds)], 0)


res = Parallel(n_jobs=a.n_jobs)(delayed(run)(e, g) for e in ENCS for g in GENES)
rows, oofs = [], []
for enc, gene, pts, y, p in res:
    m = regression(y, p)
    ci = {k: stats.paired_bootstrap(lambda i, k=k: regression(y[i], p[i])[k], len(y), a.boot) for k in m}
    rows.append({"encoder": enc, "label": registry.label(enc), "gene": gene, "n": len(y),
                 **{k: v for k, v in m.items()}, **{f"{k}_lo": ci[k][0] for k in m}, **{f"{k}_hi": ci[k][1] for k in m}})
    oofs.append(pd.DataFrame({"encoder": enc, "gene": gene, "patient": pts, "y": y, "pred": p}))
R = pd.DataFrame(rows)
R.to_csv(os.path.join(a.out, "regression.csv"), index=False)
O = pd.concat(oofs)
O.to_csv(os.path.join(a.out, "oof_reg.csv.gz"), index=False)
log("regression (Pearson r [95% CI], R²):\n" + R.assign(
    r=lambda x: x.pearson.round(3).astype(str) + " [" + x.pearson_lo.round(3).astype(str) + ", " +
    x.pearson_hi.round(3).astype(str) + "]").pivot(index="encoder", columns="gene", values="r").to_string())

# regression -> classification, vs the classifier trained on each label (X1, same encoder, same patients)
x1 = pd.read_csv(os.path.join(a.x1, "oof.csv.gz"))
crows = []
for enc, gene, pts, y, p in res:
    pr = pd.Series(p, index=pts)
    for ihc, mol in GENES[gene]:
        for lab in (ihc, mol):
            c = x1[(x1.encoder == enc) & (x1.task == lab) & x1.group.str.contains(ihc)].set_index("patient")
            ids = c.index.intersection(pr.index)
            yl = c.loc[ids, "y"].astype(int).to_numpy()
            dl = stats.delong(yl, c.loc[ids, "p"].to_numpy(), yl, pr.loc[ids].to_numpy())
            crows.append({"encoder": enc, "gene": gene, "label": lab, "n": len(ids), "auc_classifier": dl["auc_ref"],
                          "auc_regression": dl["auc_alt"], "diff": dl["diff"], "lo": dl["lo"], "hi": dl["hi"],
                          "p": dl["p"], "auc_measured_mRNA": stats.safe_auc(yl, L.loc[ids, gene].to_numpy())})
C = pd.DataFrame(crows)
C.to_csv(os.path.join(a.out, "reg_to_clf.csv"), index=False)
log("regression score vs classifier (AUC):\n" + C.round(3).to_string(index=False))
log("BENCH_REGRESSION_OK")
