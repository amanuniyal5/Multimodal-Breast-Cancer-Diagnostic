"""Week 3: label-free audit with zero-shot CONCH. Does a model that never saw an IHC or RNA label agree more with RNA?

CONCH v1 is a vision-language model: image and text embeddings live in one space, so a patch can be scored against
text prompts ("ER-positive breast carcinoma" vs "ER-negative breast carcinoma") without any training label.
We score the existing CONCH patch bags (TRACE conch_v1 20x, ≤ 4,096 random patches per patient, DX slides).

--check   Day-1 check: which CONCH space are the TRACE features in? TRIDENT/CLAM extract CONCH with
          encode_image(proj_contrast=False) (512-d, BEFORE the contrastive projection); zero-shot needs the projected
          space. (1) find the projection inside the CONCH model (checked on random images); (2) run TRIDENT's
          conch_v1 encoder on the same images to see which space it returns; (3) compare our TRIDENT re-extraction of
          a few DX slides (--ours, from extract_frozen.py --encoder conch_v1) with the TRACE patient means, raw and
          projected. Writes results/conch_feature_check.json ({"trace_space": "pre" | "post", ...}).
(default) Zero-shot audit (uses the space recorded by --check, or --space):
  scores    per patch: prompt-ensemble text embeddings (22 CONCH templates × class names), cosine × logit scale,
            softmax over the two classes; per patient: top-k pooling (mean of the k = 100 highest patch logits per
            class, MI-Zero style; PRIMARY, stated before looking) and mean of patch probabilities. Prompt sets for ER:
            all, receptor words only, subtype words only (luminal / basal-like).
  1a        AUC of the zero-shot score against each reference (no training, so no folds): paired ΔAUC
            molecular − IHC on the same patients, paired bootstrap 95% CI (2,000).
  1b        ER discordant patients (IHC ≠ RNA): side taken at the Youden cut-off chosen on CONCORDANT patients
            (binomial test, Wilson CI), and threshold-free: AUC_disc = P(score of an IHC−/RNA+ patient > score of an
            IHC+/RNA− patient); > 0.5 means the score follows RNA (Mann–Whitney p, bootstrap CI).
  concepts  14 morphology concepts (tubules, pleomorphism, mitoses, necrosis, TILs, apocrine, ...): per patch softmax
            over concepts, per patient the mean ("share of the slide that looks like the concept").
            Concept-only model: logistic regression (C = 1, balanced) on the 14 standardised concept shares, fixed
            folds, 5 seeds: 1a (IHC vs RNA label) and 1b (train concordant, score discordant).
            Concept profile by group (concordant ER+/ER−, discordant by direction and PAM50): mean z-scores.

Usage: python scripts/conch_zeroshot.py --check --ours data/features/dx_check_conch_v1
       python scripts/conch_zeroshot.py
Outputs: results/conch_zeroshot_scores.csv (per patient), conch_zeroshot_1a.csv, conch_zeroshot_1b.csv,
         conch_concepts_model.csv, conch_concepts_by_group.csv, conch_concepts_coef.csv
"""
import argparse
import glob
import json
import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.mixture import GaussianMixture
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from common import site_folds

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--feats", default="data/features")
ap.add_argument("--name", default="conch_v1_20x")
ap.add_argument("--cap", type=int, default=4096)
ap.add_argument("--labels", default="data/metadata/labels_patient.csv")
ap.add_argument("--check", action="store_true")
ap.add_argument("--ours", default="data/features/dx_check_conch_v1")
ap.add_argument("--space", default="auto", choices=["auto", "pre", "post"])
ap.add_argument("--topk", type=int, default=100)
ap.add_argument("--seeds", type=int, default=5)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--out", default="results")
a = ap.parse_args()
R = lambda f: os.path.join(a.out, f)  # noqa: E731
os.makedirs(a.out, exist_ok=True)

import torch  # noqa: E402
import torch.nn.functional as F  # noqa: E402
from conch.open_clip_custom import create_model_from_pretrained, get_tokenizer, tokenize  # noqa: E402

torch.set_grad_enabled(False)
dev = "cuda" if torch.cuda.is_available() else "cpu"
model, _ = create_model_from_pretrained("conch_ViT-B-16", "hf_hub:MahmoodLab/conch")
model = model.to(dev).eval()


def find_projection():
    """the map from encode_image(proj_contrast=False) to encode_image(proj_contrast=True), both unnormalised"""
    g = torch.Generator().manual_seed(0)
    x = torch.randn(8, 3, 448, 448, generator=g).to(dev)
    pre = model.encode_image(x, proj_contrast=False, normalize=False)
    post = model.encode_image(x, proj_contrast=True, normalize=False)
    for name, m in model.visual.named_modules():
        if name and isinstance(m, (torch.nn.Linear, torch.nn.Sequential)):
            try:
                y = m(pre)
            except Exception:  # noqa: BLE001
                continue
            if y.shape == post.shape and torch.allclose(y, post, atol=1e-3, rtol=1e-3):
                return m, f"module visual.{name}", x, pre, post
    for name, p in model.visual.named_parameters():
        if p.ndim == 2 and p.shape[0] == pre.shape[1] and p.shape[1] == post.shape[1] \
                and torch.allclose(pre @ p, post, atol=1e-3, rtol=1e-3):
            return (lambda z, p=p: z @ p), f"parameter visual.{name}", x, pre, post
    raise SystemExit(f"projection not found; model.visual = {model.visual}")


proj, proj_name, x_chk, pre_chk, post_chk = find_projection()
print(f"CONCH projection: {proj_name}; pre {tuple(pre_chk.shape)} -> post {tuple(post_chk.shape)}", flush=True)


def cos(u, v):
    return (u * v).sum(1) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1)


if a.check:
    res = {"projection": proj_name}
    # (2) which space does TRIDENT's conch_v1 encoder return?
    try:
        import inspect

        from trident.patch_encoder_models import encoder_factory
        enc = encoder_factory("conch_v1").to(dev).eval()
        print(inspect.getsource(type(enc).forward), flush=True)
        t = enc(x_chk).float()
        res["trident_vs_pre_maxabs"] = float((t - pre_chk).abs().max())
        res["trident_vs_post_maxabs"] = float((t - post_chk).abs().max())
        res["trident_space"] = "pre" if res["trident_vs_pre_maxabs"] < res["trident_vs_post_maxabs"] else "post"
    except Exception as e:  # noqa: BLE001
        res["trident_error"] = repr(e)
    # (3) our re-extraction vs TRACE patient means, raw and projected (projection is affine: mean commutes)
    z = np.load(os.path.join(a.feats, f"{a.name}_mean.npz"), allow_pickle=True)
    TRACE = pd.DataFrame(z["X"], index=z["patients"])
    files = sorted(glob.glob(os.path.join(a.ours, "*.npz")))
    ours = pd.DataFrame({os.path.basename(f)[:12]: np.load(f)["mean"] for f in files}).T
    ours = ours[ours.index.isin(TRACE.index)]
    t = TRACE.loc[ours.index].values
    o_proj = proj(torch.tensor(ours.values, dtype=torch.float32, device=dev)).cpu().numpy()
    other = TRACE.drop(ours.index).sample(len(ours), random_state=0).values
    res.update(n_slides=len(ours), cos_raw_median=float(np.median(cos(ours.values, t))),
               cos_raw_min=float(cos(ours.values, t).min()), cos_proj_median=float(np.median(cos(o_proj, t))),
               cos_raw_other_patient=float(np.median(cos(ours.values, other))))
    ours_space = res.get("trident_space", "pre")
    if res["cos_raw_median"] > 0.99:
        res["trace_space"] = ours_space
    elif res["cos_proj_median"] > 0.99 and ours_space == "pre":
        res["trace_space"] = "post"
    else:
        res["trace_space"] = "unclear"
    print(json.dumps(res, indent=1))
    json.dump(res, open(R("conch_feature_check.json"), "w"), indent=1)
    raise SystemExit(0 if res["trace_space"] in ("pre", "post") else "TRACE CONCH space unclear: do not run zero-shot")

space = a.space
if space == "auto":
    space = json.load(open(R("conch_feature_check.json")))["trace_space"]
    assert space in ("pre", "post"), space
print(f"TRACE CONCH features are in the {space}-projection space", flush=True)

# ---------------- prompts ----------------
TEMPLATES = ["CLASSNAME.", "a photomicrograph showing CLASSNAME.", "a photomicrograph of CLASSNAME.",
             "an image of CLASSNAME.", "an image showing CLASSNAME.", "an example of CLASSNAME.",
             "CLASSNAME is shown.", "this is CLASSNAME.", "there is CLASSNAME.",
             "a histopathological image showing CLASSNAME.", "a histopathological image of CLASSNAME.",
             "a histopathological photograph of CLASSNAME.", "a histopathological photograph showing CLASSNAME.",
             "shows CLASSNAME.", "presence of CLASSNAME.", "CLASSNAME is present.", "an H&E stained image of CLASSNAME.",
             "an H&E stained image showing CLASSNAME.", "an H&E image showing CLASSNAME.", "an H&E image of CLASSNAME.",
             "CLASSNAME, H&E stain.", "CLASSNAME, H&E."]
TASKS = {  # (positive class names, negative class names)
    "ER_all": (["estrogen receptor positive breast cancer", "ER-positive breast carcinoma",
                "hormone receptor positive breast cancer", "luminal breast cancer", "luminal-type invasive breast carcinoma"],
               ["estrogen receptor negative breast cancer", "ER-negative breast carcinoma",
                "hormone receptor negative breast cancer", "basal-like breast cancer", "basal-like invasive breast carcinoma"]),
    "ER_receptor": (["estrogen receptor positive breast cancer", "ER-positive breast carcinoma",
                     "hormone receptor positive breast cancer"],
                    ["estrogen receptor negative breast cancer", "ER-negative breast carcinoma",
                     "hormone receptor negative breast cancer"]),
    "ER_subtype": (["luminal breast cancer", "luminal-type invasive breast carcinoma"],
                   ["basal-like breast cancer", "basal-like invasive breast carcinoma"]),
    "PR": (["progesterone receptor positive breast cancer", "PR-positive breast carcinoma"],
           ["progesterone receptor negative breast cancer", "PR-negative breast carcinoma"]),
    "HER2": (["HER2-positive breast cancer", "HER2-amplified breast carcinoma", "HER2-enriched breast cancer"],
             ["HER2-negative breast cancer", "HER2-non-amplified breast carcinoma"]),
}
CONCEPTS = {
    "tubules": ["well-formed tubules", "tubule formation in invasive carcinoma", "glandular structures"],
    "low_grade_nuclei": ["low-grade nuclei with mild pleomorphism", "small uniform nuclei"],
    "pleomorphism": ["marked nuclear pleomorphism", "high-grade pleomorphic nuclei",
                     "large vesicular nuclei with prominent nucleoli"],
    "mitoses": ["numerous mitotic figures", "mitoses"],
    "necrosis": ["tumor necrosis", "comedo necrosis", "geographic necrosis"],
    "lymphocytes": ["tumor-infiltrating lymphocytes", "dense lymphocytic infiltrate"],
    "apocrine": ["apocrine differentiation", "apocrine cells with abundant eosinophilic granular cytoplasm"],
    "solid_sheets": ["solid sheets of tumor cells", "syncytial growth of tumor cells"],
    "lobular": ["single-file infiltration of discohesive cells", "invasive lobular carcinoma"],
    "mucin": ["extracellular mucin", "mucinous carcinoma"],
    "stroma": ["desmoplastic stroma", "fibrous stroma"],
    "adipose": ["adipose tissue", "fat cells"],
    "dcis": ["ductal carcinoma in situ"],
    "normal_breast": ["normal breast ducts and lobules", "benign breast tissue"],
}
tok = get_tokenizer()


def text_embedding(names):
    """prompt ensemble: mean of normalised embeddings over templates and class names, renormalised"""
    prompts = [t.replace("CLASSNAME", n) for n in names for t in TEMPLATES]
    e = model.encode_text(tokenize(texts=prompts, tokenizer=tok).to(dev))
    e = F.normalize(e.float(), dim=-1).mean(0)
    return F.normalize(e, dim=-1)


T_TASK = {k: torch.stack([text_embedding(p), text_embedding(n)]) for k, (p, n) in TASKS.items()}
T_CON = torch.stack([text_embedding(v) for v in CONCEPTS.values()])
scale = model.logit_scale.exp().item()
print(f"text embeddings ready ({len(TASKS)} tasks, {len(CONCEPTS)} concepts); logit scale {scale:.1f}", flush=True)

# ---------------- per-patient scores ----------------
idx = pd.read_csv(os.path.join(a.feats, f"{a.name}_index.csv")).set_index("patient")
BAGS = np.load(os.path.join(a.feats, f"{a.name}_cap{a.cap}.npy"), mmap_mode="r")
rows = []
for i, p in enumerate(idx.index):
    x = torch.tensor(np.asarray(BAGS[idx.at[p, "row"], :idx.at[p, "n_kept"]]), dtype=torch.float32, device=dev)
    if space == "pre":
        x = proj(x)
    x = F.normalize(x, dim=-1)
    r = {"patient": p, "n_patches": len(x)}
    for k, T in T_TASK.items():
        lg = scale * x @ T.T  # patches × 2 (positive, negative)
        r[f"{k}_mean"] = torch.softmax(lg, 1)[:, 0].mean().item()
        kk = min(a.topk, len(x))
        top = lg.topk(kk, dim=0).values.mean(0)
        r[f"{k}_topk"] = torch.softmax(top, 0)[0].item()
    pc = torch.softmax(scale * x @ T_CON.T, 1).mean(0).cpu().numpy()
    r.update({f"c_{c}": v for c, v in zip(CONCEPTS, pc)})
    rows.append(r)
    if i % 200 == 0:
        print(f"scored {i + 1}/{len(idx)} patients", flush=True)
S = pd.DataFrame(rows).set_index("patient")
S.to_csv(R("conch_zeroshot_scores.csv"))

# ---------------- labels ----------------
L = pd.read_csv(a.labels, index_col=0)


def gmm_pos(x):  # as in train_mil.py
    x = x.dropna()
    g = GaussianMixture(2, random_state=0).fit(x.values.reshape(-1, 1))
    return pd.Series(g.predict(x.values.reshape(-1, 1)) == np.argmax(g.means_.ravel()), index=x.index).astype(float)


L["ER_RPPApos"] = gmm_pos(L.ER_RPPA).reindex(L.index)
L["HER2_RPPApos"] = gmm_pos(L.HER2_RPPA).reindex(L.index)
q = L.PGR_log2.quantile(1 - L.PR_IHC.mean())
L["PGR_prevmatch"] = (L.PGR_log2 >= q).astype(float).where(L.PGR_log2.notna())
L = L.loc[L.index.intersection(S.index)]
rng = np.random.default_rng(0)


def paired_delta(y1, y2, s, boot):
    """AUC(y2) − AUC(y1) for one score s on the same patients; paired bootstrap CI and two-sided p"""
    est = roc_auc_score(y2, s) - roc_auc_score(y1, s)
    bs = []
    for _ in range(boot):
        i = rng.integers(0, len(s), len(s))
        if len(set(y1[i])) == 2 and len(set(y2[i])) == 2:
            bs.append(roc_auc_score(y2[i], s[i]) - roc_auc_score(y1[i], s[i]))
    bs = np.array(bs)
    return est, *np.percentile(bs, [2.5, 97.5]), 2 * min((bs <= 0).mean(), (bs >= 0).mean())


# ---------------- 1a zero-shot ----------------
COMP = [("ER", "ER_IHC", "ESR1_RNApos"), ("ER", "ER_IHC", "ER_RPPApos"), ("PR", "PR_IHC", "PGR_RNApos"),
        ("PR", "PR_IHC", "PGR_prevmatch"), ("HER2", "HER2_STRICT", "ERBB2_AMP"), ("HER2", "HER2_STRICT", "HER2_RPPApos")]
rows = []
for marker, ihc, mol in COMP:
    d = L[[ihc, mol]].dropna()
    y1, y2 = d[ihc].astype(int).values, d[mol].astype(int).values
    for task in [t for t in TASKS if t.startswith(marker)]:
        for pool in ("topk", "mean"):
            s = S.loc[d.index, f"{task}_{pool}"].values
            est, lo, hi, p = paired_delta(y1, y2, s, a.boot)
            rows.append({"marker": marker, "comparison": f"{mol} - {ihc}", "prompts": task, "pooling": pool,
                         "primary": task in ("ER_all", "PR", "HER2") and pool == "topk", "n": len(d),
                         "AUC_ihc": roc_auc_score(y1, s), "AUC_mol": roc_auc_score(y2, s), "delta": est,
                         "delta_lo": lo, "delta_hi": hi, "p_boot": p})
            print(f"1a zero-shot {task:11s} {pool:4s} {mol} - {ihc}: n={len(d)} AUC IHC {rows[-1]['AUC_ihc']:.3f} "
                  f"mol {rows[-1]['AUC_mol']:.3f} delta {est:+.3f} [{lo:+.3f}, {hi:+.3f}] p={p:.3f}", flush=True)
pd.DataFrame(rows).to_csv(R("conch_zeroshot_1a.csv"), index=False)


# ---------------- 1b zero-shot ----------------
def wilson(k, n, z=1.96):
    ph = k / n
    c = (ph + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def side(score, d, ihc, mol, label):
    """score: Series over d.index. Youden cut-off on concordant patients; threshold-free AUC_disc"""
    conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
    fpr, tpr, thr = roc_curve(conc[ihc].astype(int), score[conc.index])
    cut = thr[np.argmax(tpr - fpr)]
    k = int(((score[disc.index] >= cut).astype(int) == disc[mol].astype(int)).sum())
    up = score[disc.index[(disc[ihc] == 0) & (disc[mol] == 1)]].values  # IHC−/RNA+
    dn = score[disc.index[(disc[ihc] == 1) & (disc[mol] == 0)]].values  # IHC+/RNA−
    auc_d = np.mean([(u > v) + 0.5 * (u == v) for u in up for v in dn])
    bs = [np.mean(rng.choice(up, len(up))[:, None] > rng.choice(dn, len(dn))[None, :]) for _ in range(a.boot)]
    w = wilson(k, len(disc))
    return {"pair": f"{ihc} vs {mol}", "model": label, "n": len(disc), "n_IHCneg_RNApos": len(up),
            "n_IHCpos_RNAneg": len(dn), "cutoff_concordant_youden": cut, "side_molecular": k,
            "side_IHC": len(disc) - k, "prop_molecular": k / len(disc), "wilson_lo": w[0], "wilson_hi": w[1],
            "binom_p": binomtest(k, len(disc)).pvalue, "AUC_disc": auc_d, "AUC_disc_lo": np.percentile(bs, 2.5),
            "AUC_disc_hi": np.percentile(bs, 97.5), "mannwhitney_p": mannwhitneyu(up, dn).pvalue}


rows = []
for ihc, mol, marker in [("ER_IHC", "ESR1_RNApos", "ER"), ("HER2_STRICT", "ERBB2_AMP", "HER2"),
                         ("PR_IHC", "PGR_prevmatch", "PR")]:
    d = L[[ihc, mol]].dropna()
    for task in [t for t in TASKS if t.startswith(marker)]:
        for pool in ("topk", "mean"):
            rows.append({**side(S[f"{task}_{pool}"], d, ihc, mol, f"zero-shot {task} {pool}"),
                         "primary": task in ("ER_all", "PR", "HER2") and pool == "topk"})
            r = rows[-1]
            print(f"1b {r['model']:28s} {r['pair']}: n={r['n']} molecular {r['side_molecular']}:{r['side_IHC']} "
                  f"(p={r['binom_p']:.3g}) | AUC_disc {r['AUC_disc']:.3f} [{r['AUC_disc_lo']:.3f}, "
                  f"{r['AUC_disc_hi']:.3f}] (MW p={r['mannwhitney_p']:.3g})", flush=True)

# ---------------- concept-only model ----------------
CC = [f"c_{c}" for c in CONCEPTS]
Xc = S[CC]


def concept_model():
    return make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=3000, class_weight="balanced"))


mrows, coefs = [], []
for marker, ihc, mol in [("ER", "ER_IHC", "ESR1_RNApos"), ("ER", "ER_IHC", "ER_RPPApos"),
                         ("PR", "PR_IHC", "PGR_prevmatch"), ("HER2", "HER2_STRICT", "ERBB2_AMP")]:
    d = L[[ihc, mol, "site"]].dropna()
    pts, X = d.index.to_numpy(), Xc.loc[d.index].values
    P = {}
    for lab in (ihc, mol):
        y = d[lab].astype(int).values
        pr = np.zeros(len(pts))
        for s in range(a.seeds):
            for tr, te in site_folds(pts, s):
                pr[te] += concept_model().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1] / a.seeds
        P[lab] = pr
        m = concept_model().fit(X, y)
        coefs.append(pd.Series(m[-1].coef_[0], index=list(CONCEPTS), name=f"{lab} ({marker})"))
    y1, y2 = d[ihc].astype(int).values, d[mol].astype(int).values
    est = roc_auc_score(y2, P[mol]) - roc_auc_score(y1, P[ihc])
    bs = []
    for _ in range(a.boot):
        i = rng.integers(0, len(pts), len(pts))
        if len(set(y1[i])) == 2 and len(set(y2[i])) == 2:
            bs.append(roc_auc_score(y2[i], P[mol][i]) - roc_auc_score(y1[i], P[ihc][i]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    mrows.append({"marker": marker, "comparison": f"{mol} - {ihc}", "n": len(pts), "AUC_ihc": roc_auc_score(y1, P[ihc]),
                  "AUC_mol": roc_auc_score(y2, P[mol]), "delta": est, "delta_lo": lo, "delta_hi": hi,
                  "p_boot": 2 * min((np.array(bs) <= 0).mean(), (np.array(bs) >= 0).mean())})
    print(f"1a concept-only {mol} - {ihc}: n={len(pts)} AUC IHC {mrows[-1]['AUC_ihc']:.3f} mol "
          f"{mrows[-1]['AUC_mol']:.3f} delta {est:+.3f} [{lo:+.3f}, {hi:+.3f}]", flush=True)
# 1b with the concept-only model (train on concordant, score discordant)
for ihc, mol in [("ER_IHC", "ESR1_RNApos"), ("HER2_STRICT", "ERBB2_AMP")]:
    d = L[[ihc, mol]].dropna()
    conc, disc = d[d[ihc] == d[mol]], d[d[ihc] != d[mol]]
    p = concept_model().fit(Xc.loc[conc.index].values, conc[ihc].astype(int)).predict_proba(Xc.loc[disc.index].values)[:, 1]
    k = int(((p >= 0.5).astype(int) == disc[mol].astype(int).values).sum())
    w = wilson(k, len(disc))
    rows.append({"pair": f"{ihc} vs {mol}", "model": "concept-only LR", "n": len(disc), "side_molecular": k,
                 "side_IHC": len(disc) - k, "prop_molecular": k / len(disc), "wilson_lo": w[0], "wilson_hi": w[1],
                 "binom_p": binomtest(k, len(disc)).pvalue})
    print(f"1b concept-only {ihc} vs {mol}: n={len(disc)} molecular {k}:{len(disc) - k} "
          f"(p={rows[-1]['binom_p']:.3g})", flush=True)
pd.DataFrame(rows).to_csv(R("conch_zeroshot_1b.csv"), index=False)
pd.DataFrame(mrows).to_csv(R("conch_concepts_model.csv"), index=False)
pd.concat(coefs, axis=1).to_csv(R("conch_concepts_coef.csv"))

# ---------------- concept profile by group (ER IHC vs RNA) ----------------
Z = (Xc - Xc.mean()) / Xc.std()
d = L[["ER_IHC", "ESR1_RNApos", "PAM50"]].dropna(subset=["ER_IHC", "ESR1_RNApos"])
grp = pd.Series("other", index=d.index)
grp[(d.ER_IHC == 1) & (d.ESR1_RNApos == 1)] = "concordant ER+"
grp[(d.ER_IHC == 0) & (d.ESR1_RNApos == 0)] = "concordant ER-"
disc_up = (d.ER_IHC == 1) & (d.ESR1_RNApos == 0)
grp[disc_up] = "IHC+/RNA- " + d.PAM50[disc_up].fillna("NA")
grp[(d.ER_IHC == 0) & (d.ESR1_RNApos == 1)] = "IHC-/RNA+"
prof = Z.loc[d.index].groupby(grp).mean()
prof.insert(0, "n", grp.value_counts())
prof.to_csv(R("conch_concepts_by_group.csv"))
print(prof.round(2).to_string(), flush=True)
h, b = grp[grp == "IHC+/RNA- Her2"].index, grp[grp == "IHC+/RNA- Basal"].index
if len(h) >= 3 and len(b) >= 3:  # HER2-enriched exception: more apocrine than the Basal discordant tumours? (lead only)
    print(f"apocrine share, IHC+/RNA- HER2-enriched (n={len(h)}) vs Basal (n={len(b)}): Mann-Whitney p = "
          f"{mannwhitneyu(Xc.loc[h, 'c_apocrine'], Xc.loc[b, 'c_apocrine']).pvalue:.3g}", flush=True)
print("saved", ", ".join(R(f) for f in ("conch_zeroshot_scores.csv", "conch_zeroshot_1a.csv", "conch_zeroshot_1b.csv",
                                         "conch_concepts_model.csv", "conch_concepts_by_group.csv",
                                         "conch_concepts_coef.csv")))
