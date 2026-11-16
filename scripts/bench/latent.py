"""Latent-class model for binary tests without a gold standard (Dawid & Skene 1979; two classes, conditional
independence given the true class). Used by experiment C (evaluation against a latent truth) and X4d (soft labels,
fitted on training patients only).

dawid_skene(R): R is an (n patients x k tests) array of 0 / 1 / nan (missing allowed).
Returns prevalence, per-test sensitivity and specificity, and each patient's posterior P(truly positive).
"""
import numpy as np


def dawid_skene(R, iters=1000, tol=1e-9, eps=1e-4):
    R = np.asarray(R, dtype=float)
    obs = ~np.isnan(R)
    X = np.nan_to_num(R)
    post = np.nanmean(R, axis=1)  # start: share of positive tests per patient
    post = np.where(np.isnan(post), 0.5, post)
    for it in range(iters):
        pi = np.clip(post.mean(), eps, 1 - eps)
        se = np.clip((post[:, None] * X * obs).sum(0) / np.maximum((post[:, None] * obs).sum(0), eps), eps, 1 - eps)
        sp = np.clip(((1 - post)[:, None] * (1 - X) * obs).sum(0) / np.maximum(((1 - post)[:, None] * obs).sum(0), eps),
                     eps, 1 - eps)
        l1 = np.log(pi) + (obs * (X * np.log(se) + (1 - X) * np.log(1 - se))).sum(1)
        l0 = np.log(1 - pi) + (obs * ((1 - X) * np.log(sp) + X * np.log(1 - sp))).sum(1)
        new = 1 / (1 + np.exp(np.clip(l0 - l1, -50, 50)))
        done = np.max(np.abs(new - post)) < tol
        post = new
        if done:
            break
    ll = np.sum(np.logaddexp(l1, l0))
    return {"prevalence": float(pi), "sens": se, "spec": sp, "posterior": post, "loglik": float(ll), "iters": it + 1}


def posterior(R, fit):
    """posterior P(positive) for new patients R under a fitted model (sens / spec / prevalence held fixed)"""
    R = np.asarray(R, dtype=float)
    obs = ~np.isnan(R)
    X = np.nan_to_num(R)
    se, sp, pi = fit["sens"], fit["spec"], fit["prevalence"]
    l1 = np.log(pi) + (obs * (X * np.log(se) + (1 - X) * np.log(1 - se))).sum(1)
    l0 = np.log(1 - pi) + (obs * ((1 - X) * np.log(sp) + X * np.log(1 - sp))).sum(1)
    return 1 / (1 + np.exp(np.clip(l0 - l1, -50, 50)))
