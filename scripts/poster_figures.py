"""Poster figures (Nov 20 showcase) from the committed result files. Sizes match the template panels (inches).
Usage: python scripts/poster_figures.py [--out poster/figs]
Needs: numpy pandas matplotlib scikit-learn   (t-SNE panel also needs data/features/{uni_v1_20x,conch_v1_20x}_mean.npz
and titan_brca_20x.parquet; skipped if missing)
"""
import argparse
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--out", default="poster/figs")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
B = "results/bench"

# colour-blind-safe palette (Okabe-Ito); one meaning per colour across all figures
BLUE, ORANGE, GREEN, RED, GREY, SKY, PURPLE = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#7F7F7F", "#56B4E9", "#CC79A7"
FAMILY = {"patch": BLUE, "slide": GREEN, "resnet": RED}  # encoder families
plt.rcParams.update({"font.size": 22, "axes.titlesize": 26, "axes.labelsize": 22, "xtick.labelsize": 20,
                     "ytick.labelsize": 20, "legend.fontsize": 19, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 100,
                     "savefig.dpi": 300, "savefig.bbox": "tight", "font.family": "DejaVu Sans"})
NAME = {"uni_v1": "UNI", "uni_v2": "UNI2-h", "gigapath": "GigaPath", "hoptimus0": "H-optimus-0", "hibou_l": "Hibou-L",
        "conch_v1": "CONCH", "uni_v1_10x": "UNI 10x", "titan": "TITAN", "gigapath_slide": "GigaPath-slide",
        "ctranspath": "CTransPath †", "chief": "CHIEF †", "resnet50": "ResNet-50 (ImageNet)"}
KIND = {k: "slide" if k in ("titan", "gigapath_slide", "chief", "conch_v1") else "patch" for k in NAME}
KIND["resnet50"] = "resnet"


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(a.out, f"{name}.{ext}"))
    plt.close(fig)
    print("saved", name)


# ---- F1  Do we need pathology foundation models?  ResNet-50 vs UNI vs TITAN (panel: Our Solution, 14 x 3.6) ----
A = pd.read_csv(f"{B}/x1_subset/aucs.csv")
tasks = [("ER_IHC", "ER (IHC)"), ("ESR1_RNApos", "ER (RNA)"), ("PR_IHC", "PR (IHC)"), ("HER2_STRICT", "HER2")]
fig, ax = plt.subplots(figsize=(14, 4.6))
w = 0.26
for i, (enc, col) in enumerate([("resnet50", RED), ("uni_v1", BLUE), ("titan", GREEN)]):
    d = A[A.encoder == enc].set_index("task").loc[[t for t, _ in tasks]]
    ax.bar(np.arange(4) + (i - 1) * w, d.auc, w, color=col, label=NAME[enc],
           yerr=[d.auc - d.lo, d.hi - d.auc], capsize=4, error_kw={"lw": 2})
ax.axhline(0.5, color="k", ls="--", lw=1.5)
ax.text(3.45, 0.505, "chance", ha="right", va="bottom", fontsize=17)
ax.set_xticks(range(4), [l for _, l in tasks])
ax.set_ylim(0.4, 1.0)
ax.set_ylabel("AUC (95% CI)")
ax.set_title("ImageNet features lose 0.11–0.18 AUC: pathology models are needed", loc="left")
ax.legend(ncol=3, loc="upper right", frameon=False)
save(fig, "f1_resnet_vs_fm")

# ---- F2  Eleven encoders: ER AUC and hospital identifiability (panel: How it works / What we learned, 7.4 wide) ----
S = pd.read_csv(f"{B}/x1/site.csv").set_index("encoder")
X = pd.read_csv(f"{B}/x1/aucs.csv")
er = X[(X.task == "ER_IHC") & ~X.encoder.str.endswith("+sd")].set_index("encoder").sort_values("auc")
fig, ax = plt.subplots(1, 2, figsize=(14, 7.2), sharey=True, gridspec_kw={"wspace": 0.08})
y = np.arange(len(er))
for i, (e, r) in enumerate(er.iterrows()):
    ax[0].errorbar(r.auc, i, xerr=[[r.auc - r.lo], [r.hi - r.auc]], fmt="o", ms=11, color=FAMILY[KIND[e]], lw=2.5)
    ax[1].barh(i, 100 * S.loc[e, "accuracy"], color=FAMILY[KIND[e]], height=0.65)
ax[0].set_yticks(y, [NAME[e] for e in er.index])
ax[0].set_xlabel("ER (IHC) AUC")
ax[0].set_xlim(0.78, 0.92)
ax[1].set_xlabel("Hospital identified from features (%)")
ax[1].axvline(100 * S.chance.iloc[0], color="k", ls="--", lw=1.5)
ax[1].set_xlim(0, 100)
ax[1].text(100 * S.chance.iloc[0] + 1, -0.9, "chance", fontsize=16)
ax[0].set_title("Accuracy is the same…", loc="left", fontsize=22)
ax[1].set_title("…hospital signature is not", loc="left", fontsize=22)
ax[1].legend(handles=[plt.Line2D([], [], marker="s", ls="", color=FAMILY[k], ms=14, label=l)
                      for k, l in (("patch", "patch encoder"), ("slide", "slide / vision-language"))],
             loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2, frameon=False)
save(fig, "f2_encoders_er_vs_hospital")

# ---- F3  Label efficiency (panel: What we learned) ----
C = pd.read_csv(f"{B}/x1/curves.csv")
C = C[C.task == "ER_IHC"]
C["n"] = C.n_train.replace("all", np.nan).astype(float)
nall = 640  # ~ training patients per fold for ER (801 patients, 5 folds)
C["n"] = C.n.fillna(nall)
M = C.groupby(["encoder", "n"]).auc.mean().unstack(0)
fig, ax = plt.subplots(figsize=(7.2, 6.4))
for e in M.columns:
    hi = e in ("titan", "uni_v1", "conch_v1")
    ax.plot(M.index, M[e], color=FAMILY[KIND[e]] if hi else "#BBBBBB", lw=4 if hi else 1.6, marker="o" if hi else None,
            ms=8, zorder=3 if hi else 1)
    if hi:
        ax.text(M.index[-1] + 15, M[e].iloc[-1] + (0.008 if e == "titan" else -0.016 if e == "uni_v1" else -0.014),
                NAME[e], color=FAMILY[KIND[e]], fontsize=19, va="center")
ax.set_xscale("log")
ax.set_xticks([25, 50, 100, 200, 400, nall], ["25", "50", "100", "200", "400", "all"])
ax.set_xlim(20, 2200)
ax.set_xlabel("Training patients")
ax.set_ylabel("ER (IHC) AUC")
ax.set_title("TITAN with 25 patients ≈ UNI with ~150", loc="left", fontsize=20)
save(fig, "f3_label_efficiency")

# ---- F4  Improvements: mean+SD pooling and mRNA-regression target (panel: What we learned / Our Approach) ----
O = X.set_index(["encoder", "task"]).auc
sd_pairs = [("gigapath", "GigaPath"), ("uni_v1_10x", "UNI 10x"), ("hoptimus0", "H-optimus-0"), ("ctranspath", "CTransPath †"),
            ("uni_v2", "UNI2-h")]
R = pd.read_csv(f"{B}/x3/reg_to_clf.csv")
R = R[R.label == "ER_IHC"].set_index("encoder")
fig, ax = plt.subplots(1, 2, figsize=(14, 5.6))
for k, (e, l) in enumerate(sd_pairs):
    ax[0].plot([0, 1], [O[(e, "ER_IHC")], O[(e + "+sd", "ER_IHC")]], marker="o", ms=11, lw=3, color=BLUE)
    ax[0].text(1.04, O[(e + "+sd", "ER_IHC")], l, va="center", fontsize=18)
ax[0].set_xticks([0, 1], ["mean pool", "mean + SD"])
ax[0].set_xlim(-0.2, 1.8)
ax[0].set_ylabel("ER (IHC) AUC")
ax[0].set_title("A cheap fix: add the SD", loc="left", fontsize=22)
rr = R.loc[["hoptimus0", "chief", "gigapath", "gigapath_slide", "uni_v1"]]
for i, (e, r) in enumerate(rr.iterrows()):
    ax[1].plot([0, 1], [r.auc_classifier, r.auc_regression], marker="o", ms=11, lw=3, color=GREEN)
    ax[1].text(1.04, r.auc_regression, NAME[e], va="center", fontsize=18)
ax[1].set_xticks([0, 1], ["train on IHC", "train on\nESR1 mRNA"])
ax[1].set_xlim(-0.2, 1.9)
ax[1].set_title("Train on expression, test on IHC", loc="left", fontsize=22)
save(fig, "f4_improvements")

# ---- F5  What the representations contain: CKA heat map + kNN purity (panel: How it works) ----
K = pd.read_csv(f"{B}/x5/cka.csv")
enc = [e for e in NAME if e in set(K.enc1) | set(K.enc2) and e != "resnet50"]
Mx = pd.DataFrame(np.nan, index=enc, columns=enc)
for r in K.itertuples():
    Mx.loc[r.enc1, r.enc2] = Mx.loc[r.enc2, r.enc1] = r.cka
order = Mx.mean().sort_values(ascending=False).index
Mx = Mx.loc[order, order]
fig, ax = plt.subplots(figsize=(7.2, 6.8))
im = ax.imshow(Mx.values, vmin=0.3, vmax=1, cmap="viridis")
ax.set_xticks(range(len(order)), [NAME[e] for e in order], rotation=60, ha="right", fontsize=15)
ax.set_yticks(range(len(order)), [NAME[e] for e in order], fontsize=15)
for i in range(len(order)):
    for j in range(len(order)):
        v = Mx.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=11, color="w" if v < 0.65 else "k")
ax.grid(False)
fig.colorbar(im, ax=ax, fraction=0.046, label="linear CKA")
ax.set_title("Patch encoders learn the same thing;\nTITAN differs", loc="left", fontsize=20)
save(fig, "f5_cka")

# ---- F6  Hospital vs biology: kNN purity by encoder (panel: What we learned) ----
E = pd.read_csv(f"{B}/x5/embedding_metrics.csv").set_index("encoder").sort_values("knn_purity_site")
fig, ax = plt.subplots(figsize=(7.2, 6.4))
for e, r in E.iterrows():
    ax.scatter(r.knn_purity_site, r.knn_purity_ESR1_RNApos, s=260, color=FAMILY[KIND[e]], zorder=3)
    ax.annotate(NAME[e], (r.knn_purity_site, r.knn_purity_ESR1_RNApos), xytext=(7, 6), textcoords="offset points",
                fontsize=15)
ax.set_xlabel("Neighbours from the same hospital")
ax.set_ylabel("Neighbours with the same ER status")
ax.set_title("Less hospital, more biology", loc="left", fontsize=22)
save(fig, "f6_purity")

# ---- F7  Calibration: expected calibration error before / after Platt scaling (panel: What we learned) ----
Kc = pd.read_csv(f"{B}/x4/calibration.csv")
Kc = Kc[(Kc.encoder == "uni_v1") & Kc.task.isin(["ER_IHC", "PR_IHC", "HER2_STRICT", "ERBB2_AMP"])]
lab = {"ER_IHC": "ER", "PR_IHC": "PR", "HER2_STRICT": "HER2", "ERBB2_AMP": "ERBB2\namp."}
fig, ax = plt.subplots(figsize=(7.2, 5.2))
for i, (v, col) in enumerate([("raw", GREY), ("platt", BLUE)]):
    d = Kc[Kc.variant == v].set_index("task").loc[list(lab)]
    ax.bar(np.arange(4) + (i - 0.5) * 0.38, d.ece, 0.38, color=col, label="raw" if v == "raw" else "after recalibration")
ax.set_xticks(range(4), list(lab.values()))
ax.set_ylabel("Calibration error (ECE)")
ax.legend(frameon=False)
ax.set_title("Raw scores are over-confident", loc="left", fontsize=22)
save(fig, "f7_calibration")


# ---- F8  Pipeline / workflow diagram (panel: How it works, ~14 x 3.6 or 7.4 wide stacked) ----
def box(ax, x, y, w, h, text, fc, tc="w", fs=19):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec="none"))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=tc, fontsize=fs, fontweight="bold", linespacing=1.0)


fig, ax = plt.subplots(figsize=(14, 3.8))
ax.set_xlim(0, 14)
ax.set_ylim(0, 3.8)
ax.axis("off")
steps = [("H&E slide\n(FFPE or frozen)", GREY), ("Tile into\n256 px patches", GREY), ("Foundation model\n(UNI, TITAN, CONCH…)", BLUE),
         ("Aggregate patches\n(mean+SD, ABMIL…)", GREEN), ("Calibrated\nER / PR / HER2", ORANGE)]
for i, (t, c) in enumerate(steps):
    box(ax, 0.1 + i * 2.8, 1.5, 2.3, 1.5, t, c)
    if i < len(steps) - 1:
        ax.add_patch(FancyArrowPatch((2.45 + i * 2.8, 2.25), (2.85 + i * 2.8, 2.25), arrowstyle="-|>", mutation_scale=28,
                                     color="k", lw=2.5))
ax.text(7, 0.75, "Second reader: flags cases where the H&E disagrees with the IHC report or is uncertain, and shows "
        "which regions drove the call", ha="center", va="center", fontsize=19)
save(fig, "f8_pipeline")

# ---- F8v  Vertical pipeline for the narrow panels (Our Approach / How it works, 7.2 x 6.0) ----
fig, ax = plt.subplots(figsize=(7.2, 6.0))
ax.set_xlim(0, 7.2)
ax.set_ylim(0, 6.0)
ax.axis("off")
for i, (t, c) in enumerate(steps):
    y0 = 5.1 - i * 1.2
    box(ax, 0.2, y0, 6.8, 0.85, t, c, fs=17)
    if i < len(steps) - 1:
        ax.add_patch(FancyArrowPatch((3.6, y0 - 0.03), (3.6, y0 - 0.33), arrowstyle="-|>", mutation_scale=26, color="k", lw=2.5))
save(fig, "f8v_pipeline_vertical")

# ---- F9  X2: aggregation methods beyond ABMIL (only if the Ada runs are in results/bench/x2) ----
p = f"{B}/x2/aucs.csv"
if os.path.exists(p):
    X2 = pd.read_csv(p)
    ms = [m for m in ["mean", "meanstd", "knn25", "abmil", "clam", "transmil", "dsmil", "acmil"] if m in set(X2.method)]
    tt = [("ER_IHC", "ER (IHC)"), ("ESR1_RNApos", "ER (RNA)"), ("PR_IHC", "PR (IHC)"), ("HER2_STRICT", "HER2")]
    fig, ax = plt.subplots(figsize=(14, 5.2))
    for i, (t, l) in enumerate(tt):
        d = X2[X2.task == t].set_index("method").loc[ms]
        ax.bar(np.arange(len(ms)) + (i - 1.5) * 0.2, d.auc, 0.2, label=l, color=[BLUE, SKY, ORANGE, GREEN][i])
    ax.set_xticks(range(len(ms)), ms)
    ax.set_ylim(0.5, 1.0)
    ax.set_ylabel("AUC")
    ax.legend(ncol=4, frameon=False)
    ax.set_title("Complex MIL heads ≈ simple mean + SD pooling (UNI)", loc="left")
    save(fig, "f9_aggregation")
else:
    print("skip F9: results/bench/x2/aucs.csv not found (run scripts/analyse_x2.py after the Ada jobs)")

# ---- F10  t-SNE of patient features: UNI vs TITAN, coloured by hospital and by ER (RNA) (panel: How it works) ----
try:
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler

    from bench import registry, tasks as T

    L = T.labels()
    uni, titan = registry.load_means("uni_v1"), registry.load_means("titan")
    pts = sorted(set(uni.index) & set(titan.index) & set(L.index))
    top = L.loc[pts, "site"].value_counts().index[:8]
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(2, 2, figsize=(14, 12.5))
    for r, (e, df) in enumerate([("uni_v1", uni), ("titan", titan)]):
        Z = StandardScaler().fit_transform(df.loc[pts].values)
        Y = TSNE(perplexity=30, init="pca", random_state=0).fit_transform(PCA(50, random_state=0).fit_transform(Z))
        site = L.loc[pts, "site"].to_numpy()
        ax[r, 0].scatter(Y[:, 0], Y[:, 1], s=14, c="#DDDDDD")
        for j, s in enumerate(top):
            m = site == s
            ax[r, 0].scatter(Y[m, 0], Y[m, 1], s=22, color=cmap(j))
        v = L.loc[pts, "ESR1_RNApos"].to_numpy()
        for val, col, nm in ((1.0, BLUE, "ER+ (RNA)"), (0.0, ORANGE, "ER− (RNA)")):
            m = v == val
            ax[r, 1].scatter(Y[m, 0], Y[m, 1], s=22, color=col, label=nm, alpha=0.8)
        ax[r, 0].set_ylabel(NAME[e])
        for c in range(2):
            ax[r, c].set_xticks([])
            ax[r, c].set_yticks([])
            ax[r, c].grid(False)
    ax[0, 0].set_title("coloured by hospital", loc="left", fontsize=22)
    ax[0, 1].set_title("coloured by ER status", loc="left", fontsize=22)
    ax[0, 1].legend(markerscale=2, frameon=False, loc="upper right")
    save(fig, "f10_tsne_uni_vs_titan")
except Exception as ex:  # features not on this machine
    print("skip F10:", repr(ex))
