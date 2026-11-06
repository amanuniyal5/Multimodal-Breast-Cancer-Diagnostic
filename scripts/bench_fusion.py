"""X4a encoder fusion / ensembling and X4f calibration, mean-pool heads (docs/DESIGN_SPEC.md §5.4).

On the same label-group cohorts as X1 (so every comparison is on identical patients):
  early fusion  concatenated per-encoder means (each z-scored inside the pipeline), logistic regression with C chosen
                on inner folds, 5 seeds
  late fusion   mean logit of the single-encoder X1 out-of-fold probabilities
  compared (DeLong) with the reference encoder and with the best single encoder (best chosen on the same data: an
  optimistic comparator, so a fusion gain over it is conservative)
Calibration (seed 0): raw class-balanced linear-probe probabilities vs Platt recalibration fitted on inner-fold
predictions of the training patients only; Brier, ECE (15 bins), calibration slope / intercept.
Outputs fusion.csv, calibration.csv, oof_fusion.csv.gz.
Usage: python scripts/bench_fusion.py --x1 results/bench/x1 --out results/bench/x4 [--fuse uni_v1,conch_v1,titan]
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import folds, registry, stats, tasks  # noqa: E402
from bench.metrics import classification  # noqa: E402
from common import linear_probe  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--x1", default="results/bench/x1")
ap.add_argument("--fuse", default="", help="encoders to fuse; default: all non-TCGA-pretrained encoders in X1")
ap.add_argument("--ref", default="uni_v1")
ap.add_argument("--calibrate", default="uni_v1,titan")
ap.add_argument("--cs", default="0.001,0.003,0.01,0.03")
ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--n-jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
ap.add_argument("--out", default="results/bench/x4")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()
CS = [float(c) for c in a.cs.split(",")]
TASKS = ["ER_IHC", "ESR1_RNApos", "ER_RPPApos", "PR_IHC", "PGR_prevmatch", "HER2_STRICT", "ERBB2_AMP"]


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


x1 = pd.read_csv(os.path.join(a.x1, "oof.csv.gz"))
x1 = x1[x1.task.isin(TASKS) & ~x1.encoder.str.endswith("+sd")]
FUSE = a.fuse.split(",") if a.fuse else [e for e in x1.encoder.unique() if not registry.ENC[e]["tcga_pretrained"]]
FE = {e: registry.load_means(e) for e in FUSE}
log(f"fusing {FUSE}")


def lr(c):
    return make_pipeline(StandardScaler(), LogisticRegression(C=c, max_iter=3000, class_weight="balanced"))


def early(task, group, seed):
    d = x1[(x1.task == task) & (x1.group == group) & (x1.encoder == FUSE[0])]
    pts, y = d.patient.to_numpy(), d.y.astype(int).to_numpy()
    X = np.hstack([FE[e].loc[pts].values for e in FUSE])
    p = np.zeros(len(y))
    for tr, te in folds.outer(pts, seed):
        sc = []
        for c in CS:
            q = np.zeros(len(tr))
            for itr, iva in folds.inner(pts[tr], seed):
                q[iva] = lr(c).fit(X[tr][itr], y[tr][itr]).predict_proba(X[tr][iva])[:, 1]
            sc.append(stats.safe_auc(y[tr], q))
        p[te] = lr(CS[int(np.nanargmax(sc))]).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return task, group, seed, pts, y, p


groups = x1[["task", "group"]].drop_duplicates().itertuples(index=False)
groups = [(t, g) for t, g in groups]
res = Parallel(n_jobs=a.n_jobs)(delayed(early)(t, g, s) for t, g in groups for s in range(a.seeds))
rows, oofs = [], []
logit = lambda p: np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))  # noqa: E731
for t, g in groups:
    r = [x for x in res if x[0] == t and x[1] == g]
    pts, y = r[0][3], r[0][4]
    pe = np.mean([x[5] for x in r], 0)
    single = {e: x1[(x1.task == t) & (x1.group == g) & (x1.encoder == e)].set_index("patient").loc[pts, "p"].to_numpy()
              for e in x1.encoder.unique()}
    pl = 1 / (1 + np.exp(-np.mean([logit(single[e]) for e in FUSE], 0)))
    best = max((e for e in single if not registry.ENC[e]["tcga_pretrained"]), key=lambda e: stats.safe_auc(y, single[e]))
    for name, p in (("early fusion", pe), ("late fusion", pl)):
        auc, lo, hi = stats.auc_ci(y, p)
        row = {"task": t, "group": g, "method": name, "n": len(y), "auc": auc, "lo": lo, "hi": hi}
        for tag, e in (("ref", a.ref), ("best", best)):
            dl = stats.delong(y, single[e], y, p)
            row.update({f"{tag}": e, f"auc_{tag}": dl["auc_ref"], f"diff_{tag}": dl["diff"], f"lo_{tag}": dl["lo"],
                        f"hi_{tag}": dl["hi"], f"p_{tag}": dl["p"]})
        rows.append(row)
        oofs.append(pd.DataFrame({"method": name, "task": t, "group": g, "patient": pts, "y": y, "p": p}))
F = pd.DataFrame(rows)
F.to_csv(os.path.join(a.out, "fusion.csv"), index=False)
pd.concat(oofs).to_csv(os.path.join(a.out, "oof_fusion.csv.gz"), index=False)
log("fusion:\n" + F[["task", "method", "auc", "best", "auc_best", "diff_best", "lo_best", "hi_best", "diff_ref",
                     "lo_ref", "hi_ref"]].round(3).to_string(index=False))


# ---- X4f calibration (seed 0) ----
def calib(enc, t, g):
    d = x1[(x1.task == t) & (x1.group == g) & (x1.encoder == enc)]
    pts, y = d.patient.to_numpy(), d.y.astype(int).to_numpy()
    X = registry.load_means(enc).loc[pts].values
    raw, cal = np.zeros(len(y)), np.zeros(len(y))
    for tr, te in folds.outer(pts, 0):
        q = np.zeros(len(tr))
        for itr, iva in folds.inner(pts[tr], 0):
            q[iva] = linear_probe().fit(X[tr][itr], y[tr][itr]).predict_proba(X[tr][iva])[:, 1]
        platt = LogisticRegression(C=1e6).fit(logit(q)[:, None], y[tr])
        raw[te] = linear_probe().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        cal[te] = platt.predict_proba(logit(raw[te])[:, None])[:, 1]
    out = []
    for name, p in (("raw", raw), ("platt", cal)):
        s = LogisticRegression(C=1e6).fit(logit(p)[:, None], y)
        out.append({"encoder": enc, "task": t, "variant": name, "n": len(y), "prevalence": y.mean(),
                    **classification(y, p), "slope": s.coef_[0, 0], "intercept": s.intercept_[0]})
    return out


cal = Parallel(n_jobs=a.n_jobs)(delayed(calib)(e, t, g) for e in a.calibrate.split(",") for t, g in groups)
K = pd.DataFrame([r for x in cal for r in x])
K.to_csv(os.path.join(a.out, "calibration.csv"), index=False)
log("calibration:\n" + K.round(3).to_string(index=False))
log("BENCH_FUSION_OK")
