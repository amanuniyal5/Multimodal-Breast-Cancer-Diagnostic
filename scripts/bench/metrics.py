"""Classification and regression scores."""
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import average_precision_score, brier_score_loss, r2_score, roc_auc_score


def ece(y, p, bins=15):
    y, p = np.asarray(y), np.asarray(p)
    edges = np.linspace(0, 1, bins + 1)
    b = np.clip(np.digitize(p, edges[1:-1]), 0, bins - 1)
    return sum(abs(y[b == k].mean() - p[b == k].mean()) * (b == k).mean() for k in range(bins) if (b == k).any())


def classification(y, p):
    y, p = np.asarray(y).astype(int), np.asarray(p, dtype=float)
    return {"auc": roc_auc_score(y, p), "auprc": average_precision_score(y, p), "brier": brier_score_loss(y, p),
            "ece": ece(y, p)}


def regression(y, yhat):
    y, yhat = np.asarray(y, dtype=float), np.asarray(yhat, dtype=float)
    return {"pearson": pearsonr(y, yhat)[0], "spearman": spearmanr(y, yhat)[0], "r2": r2_score(y, yhat)}
