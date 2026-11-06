"""X5a / X5c / X5e: what the encoders' slide representations contain (docs/DESIGN_SPEC.md §5.5).

X5a  t-SNE (PCA 50 -> t-SNE, perplexity 30; UMAP too if umap-learn is installed) of patient mean features per encoder,
     coloured by ER (RNA and IHC), PAM50 and hospital; the ER-discordant patients drawn as crosses. Quantities in the
     full feature space (z-scored, cosine, leave-one-out 10 nearest neighbours): label purity for ER_IHC, ESR1_RNApos,
     PAM50 and hospital; silhouette by PAM50; for each ER-discordant patient, the share of its 10 nearest CONCORDANT
     patients whose ER status equals the patient's RNA status (0.5 = no side; Wilcoxon vs 0.5).
X5c  linear CKA between encoders (same patients), bootstrap 95% CI.
X5e  concept probing: ridge from each encoder's means to the 14 CONCH concept shares (results/conch_zeroshot_scores.csv),
     out-of-fold R² on the fixed seed-0 folds.
Outputs (--out): embedding_metrics.csv, discordant_neighbours.csv, cka.csv, concept_probe.csv, figures/*.png
Usage: python scripts/embeddings.py --out results/bench/x5 [--encoders ...]
"""
import argparse
import os
import sys
import time
import warnings

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.manifold import TSNE
from sklearn.metrics import r2_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench import folds, registry, tasks  # noqa: E402

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--encoders", default="")
ap.add_argument("--k", type=int, default=10)
ap.add_argument("--boot", type=int, default=200)
ap.add_argument("--concepts", default="results/conch_zeroshot_scores.csv")
ap.add_argument("--out", default="results/bench/x5")
a = ap.parse_args()
os.makedirs(os.path.join(a.out, "figures"), exist_ok=True)
t0 = time.time()


def log(*x):
    print(f"[{time.time() - t0:6.0f}s]", *x, flush=True)


L = tasks.labels()
ENCS = a.encoders.split(",") if a.encoders else [e for e in registry.available() if e != "resnet50"]
FE = {e: registry.load_means(e) for e in ENCS}
pts = set(L.index) & set(pd.read_csv("data/splits/folds.csv").patient)
for e in ENCS:
    pts &= set(FE[e].index)
pts = sorted(pts)
L = L.loc[pts]
er = L[["ER_IHC", "ESR1_RNApos"]].dropna()
disc = er.index[er.ER_IHC != er.ESR1_RNApos]
conc = er.index[er.ER_IHC == er.ESR1_RNApos]
log(f"encoders {ENCS}; {len(pts)} patients; ER-discordant {len(disc)}")


def zs(enc):
    return StandardScaler().fit_transform(FE[enc].loc[pts].values)


# ---------------- X5a ----------------
mrows, nrows, emb2d = [], [], {}
for e in ENCS:
    Z = zs(e)
    nn = NearestNeighbors(n_neighbors=a.k + 1, metric="cosine").fit(Z)
    nb = nn.kneighbors(Z, return_distance=False)[:, 1:]
    row = {"encoder": e, "label": registry.label(e), "n": len(pts)}
    for col in ["ER_IHC", "ESR1_RNApos", "PAM50", "site"]:
        v = L[col].to_numpy()
        ok = pd.notna(v)
        same = np.array([np.mean([v[j] == v[i] for j in nb[i] if pd.notna(v[j])]) if ok[i] else np.nan
                         for i in range(len(v))])
        row[f"knn_purity_{col}"] = np.nanmean(same)
        row[f"chance_{col}"] = (pd.Series(v[ok]).value_counts(normalize=True) ** 2).sum()
    lab = L.PAM50.to_numpy()
    okp = pd.notna(lab)
    row["silhouette_PAM50"] = silhouette_score(Z[okp], lab[okp], metric="cosine")
    # discordant patients: side of their nearest concordant neighbours
    ci = np.array([pts.index(p) for p in conc])
    nnc = NearestNeighbors(n_neighbors=a.k, metric="cosine").fit(Z[ci])
    di = np.array([pts.index(p) for p in disc])
    nbc = nnc.kneighbors(Z[di], return_distance=False)
    cstat = er.loc[conc, "ER_IHC"].to_numpy()
    share = np.array([(cstat[nbc[r]] == er.loc[p, "ESR1_RNApos"]).mean() for r, p in enumerate(disc)])
    for r, p in enumerate(disc):
        nrows.append({"encoder": e, "patient": p, "ER_IHC": er.loc[p, "ER_IHC"], "ESR1_RNApos": er.loc[p, "ESR1_RNApos"],
                      "PAM50": L.loc[p, "PAM50"], "share_rna_side": share[r]})
    row.update(disc_n=len(disc), disc_share_rna_side=share.mean(), disc_majority_rna=int((share > 0.5).sum()),
               disc_majority_ihc=int((share < 0.5).sum()), disc_wilcoxon_p=wilcoxon(share - 0.5).pvalue)
    mrows.append(row)
    emb2d[e] = TSNE(perplexity=30, init="pca", random_state=0).fit_transform(PCA(50, random_state=0).fit_transform(Z))
    log(f"{e}: purity ER-RNA {row['knn_purity_ESR1_RNApos']:.3f}, ER-IHC {row['knn_purity_ER_IHC']:.3f}, "
        f"PAM50 {row['knn_purity_PAM50']:.3f}, hospital {row['knn_purity_site']:.3f}; discordant RNA-side share "
        f"{share.mean():.3f} (p {row['disc_wilcoxon_p']:.3g})")
M = pd.DataFrame(mrows)
M.to_csv(os.path.join(a.out, "embedding_metrics.csv"), index=False)
pd.DataFrame(nrows).to_csv(os.path.join(a.out, "discordant_neighbours.csv"), index=False)

try:
    import umap  # optional
    for e in ENCS:
        emb2d[e + "_umap"] = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=0).fit_transform(zs(e))
except ImportError:
    log("umap-learn not installed: t-SNE only")

PAM = {"LumA": "#3b6fb6", "LumB": "#7fb3e6", "Her2": "#d9822b", "Basal": "#c0392b", "Normal": "#7f8c8d"}
top = L.site.value_counts().index[:10]
cmap = plt.get_cmap("tab10")
dmask = np.isin(pts, disc)
for key, Y in emb2d.items():
    fig, ax = plt.subplots(1, 4, figsize=(22, 5.6))
    for i, (col, title) in enumerate([("ESR1_RNApos", "ER by RNA"), ("ER_IHC", "ER by IHC")]):
        v = L[col].to_numpy()
        for val, c, nm in ((1.0, "#2e86c1", "positive"), (0.0, "#e67e22", "negative")):
            m = (v == val) & ~dmask
            ax[i].scatter(Y[m, 0], Y[m, 1], s=5, c=c, alpha=0.5, label=nm)
            m = (v == val) & dmask
            ax[i].scatter(Y[m, 0], Y[m, 1], s=40, c=c, marker="x", linewidths=1.6)
        ax[i].set_title(f"{title} (x = IHC/RNA discordant)")
        ax[i].legend(markerscale=3, fontsize=8)
    for k, c in PAM.items():
        m = (L.PAM50 == k).to_numpy()
        ax[2].scatter(Y[m, 0], Y[m, 1], s=5, c=c, alpha=0.6, label=k)
    ax[2].scatter(Y[dmask, 0], Y[dmask, 1], s=40, c="k", marker="x", linewidths=1.2, label="discordant")
    ax[2].set_title("PAM50 subtype")
    ax[2].legend(markerscale=3, fontsize=8)
    ax[3].scatter(Y[:, 0], Y[:, 1], s=4, c="#d5d8dc")
    for j, s in enumerate(top):
        m = (L.site == s).to_numpy()
        ax[3].scatter(Y[m, 0], Y[m, 1], s=6, color=cmap(j), label=s)
    ax[3].set_title("Hospital (10 largest)")
    ax[3].legend(markerscale=3, fontsize=7, ncol=2)
    for x in ax:
        x.set_xticks([]), x.set_yticks([])
    e = key.replace("_umap", "")
    fig.suptitle(f"{registry.label(e)}: {'UMAP' if key.endswith('_umap') else 't-SNE'} of patient mean features "
                 f"(n = {len(pts)})", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "figures", f"x5a_{key}.png"), dpi=110)
    plt.close(fig)
log("embedding metrics:\n" + M.drop(columns=["label"]).round(3).to_string(index=False))


# ---------------- X5c CKA ----------------
def cka(X, Y):
    X, Y = X - X.mean(0), Y - Y.mean(0)
    return np.linalg.norm(Y.T @ X) ** 2 / (np.linalg.norm(X.T @ X) * np.linalg.norm(Y.T @ Y))


Zs = {e: zs(e) for e in ENCS}
r = np.random.default_rng(0)
boots = [r.integers(0, len(pts), len(pts)) for _ in range(a.boot)]
crows = []
for i, e1 in enumerate(ENCS):
    for e2 in ENCS[i:]:
        v = cka(Zs[e1], Zs[e2])
        b = [cka(Zs[e1][ix], Zs[e2][ix]) for ix in boots] if e1 != e2 else [1.0]
        crows.append({"enc1": e1, "enc2": e2, "cka": v, "lo": np.percentile(b, 2.5), "hi": np.percentile(b, 97.5)})
C = pd.DataFrame(crows)
C.to_csv(os.path.join(a.out, "cka.csv"), index=False)
mat = C.pivot(index="enc1", columns="enc2", values="cka").reindex(index=ENCS, columns=ENCS)
mat = mat.combine_first(mat.T)
fig, ax = plt.subplots(figsize=(1.0 * len(ENCS) + 2, 0.9 * len(ENCS) + 1.5))
im = ax.imshow(mat.values, vmin=0, vmax=1, cmap="viridis")
ax.set_xticks(range(len(ENCS)), [registry.label(e) for e in ENCS], rotation=45, ha="right")
ax.set_yticks(range(len(ENCS)), [registry.label(e) for e in ENCS])
for i in range(len(ENCS)):
    for j in range(len(ENCS)):
        ax.text(j, i, f"{mat.values[i, j]:.2f}", ha="center", va="center", color="w" if mat.values[i, j] < 0.6 else "k",
                fontsize=8)
fig.colorbar(im, ax=ax, label="linear CKA")
ax.set_title("Representation similarity between encoders (patient means)")
fig.tight_layout()
fig.savefig(os.path.join(a.out, "figures", "x5c_cka.png"), dpi=120)
plt.close(fig)
log("CKA:\n" + mat.round(2).to_string())

# ---------------- X5e concept probing ----------------
cs = pd.read_csv(a.concepts).set_index("patient")
ccols = [c for c in cs.columns if c.startswith("c_")]
ids = [p for p in pts if p in cs.index]
prow = []
for e in ENCS:
    X = FE[e].loc[ids].values
    for c in ccols:
        y = cs.loc[ids, c].to_numpy()
        pr = np.zeros(len(y))
        for tr, te in folds.outer(ids, 0):
            m = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(0, 5, 11))).fit(X[tr], y[tr])
            pr[te] = m.predict(X[te])
        prow.append({"encoder": e, "concept": c[2:], "n": len(ids), "r2": r2_score(y, pr)})
P = pd.DataFrame(prow)
P.to_csv(os.path.join(a.out, "concept_probe.csv"), index=False)
log("concept probing R²:\n" + P.pivot(index="concept", columns="encoder", values="r2").round(2).to_string())
log("EMBEDDINGS_OK")
