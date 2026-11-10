"""MIL aggregation heads for X2 (docs/DESIGN_SPEC.md §5.2): plain PyTorch, no external MIL packages.

Shared interface: forward(x) -> (logit, aux) for one bag x (n_patches, d); objective(logit, aux, y, lossf) -> loss.
aux carries attention / instance scores for the feature maps (X5b). Binary task, one logit per bag.

  abmil     gated attention MIL (Ilse et al. 2018), CLAM-style sizes; identical to the class in train_mil.py (same layer
            order, so the same seed gives the same initialisation)
  clam      CLAM-SB (Lu et al. 2021): gated attention + instance classifiers on the top-8 / bottom-8 attended patches;
            loss 0.7 bag + 0.3 instance (cross-entropy instead of the smooth top-1 SVM, a CLAM option)
  transmil  TransMIL (Shao et al. 2021): 1024 -> 512, CLS token, 2 self-attention layers (8 heads) with the PPEG
            positional convolution between them; exact attention (PyTorch memory-efficient kernel) instead of the
            Nystrom approximation
  dsmil     DSMIL (Li et al. 2021): instance classifier -> critical instance -> attention by query similarity to it;
            loss 0.5 bag + 0.5 max-instance; prediction = mean of the two logits
  acmil     ACMIL (Zhang et al. 2023): 5 attention branches with stochastic top-10 instance masking (p 0.6) during
            training, one classifier per branch, bag classifier on the branch-averaged attention, diversity loss =
            mean pairwise cosine similarity of the branch attention maps (weight 0.1)
"""
import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class MIL(nn.Module):
    def objective(self, logit, aux, y, lossf):
        return lossf(logit, y)


class ABMIL(MIL):  # gated attention MIL (Ilse et al. 2018), CLAM-style sizes: copy of train_mil.py
    def __init__(self, d, h=256, att=128, p=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
        self.a = nn.Sequential(nn.Linear(h, att), nn.Tanh())
        self.b = nn.Sequential(nn.Linear(h, att), nn.Sigmoid())
        self.c = nn.Linear(att, 1)
        self.cls = nn.Linear(h, 1)

    def forward(self, x):
        hdn = self.fc(x)
        w = torch.softmax(self.c(self.a(hdn) * self.b(hdn)).squeeze(-1), dim=0)
        return self.cls((w[:, None] * hdn).sum(0)).squeeze(), {"attention": w}


class CLAMSB(MIL):
    def __init__(self, d, h=512, att=256, p=0.25, k=8, w_bag=0.7):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
        self.a = nn.Sequential(nn.Linear(h, att), nn.Tanh(), nn.Dropout(p))
        self.b = nn.Sequential(nn.Linear(h, att), nn.Sigmoid(), nn.Dropout(p))
        self.c = nn.Linear(att, 1)
        self.cls = nn.Linear(h, 1)
        self.inst = nn.ModuleList([nn.Linear(h, 2) for _ in range(2)])  # one instance classifier per class
        self.k, self.w_bag = k, w_bag

    def forward(self, x):
        hdn = self.fc(x)
        s = self.c(self.a(hdn) * self.b(hdn)).squeeze(-1)
        w = torch.softmax(s, dim=0)
        return self.cls((w[:, None] * hdn).sum(0)).squeeze(), {"attention": w, "h": hdn, "score": s}

    def objective(self, logit, aux, y, lossf):
        h, s = aux["h"], aux["score"]
        k = max(1, min(self.k, len(s) // 2))
        top, bot = torch.topk(s, k).indices, torch.topk(-s, k).indices
        c = int(y.item())
        one, zero = torch.ones(k, dtype=torch.long, device=s.device), torch.zeros(k, dtype=torch.long, device=s.device)
        # in-class: most attended -> positive evidence, least attended -> negative; out-of-class: most attended -> negative
        l_in = F.cross_entropy(self.inst[c](torch.cat([h[top], h[bot]])).float(), torch.cat([one, zero]))
        l_out = F.cross_entropy(self.inst[1 - c](h[top]).float(), zero)
        return self.w_bag * lossf(logit, y) + (1 - self.w_bag) * 0.5 * (l_in + l_out)


class _SelfAttention(nn.Module):
    def __init__(self, dim, heads=8, p=0.1):
        super().__init__()
        self.norm = nn.LayerNorm(dim)
        self.qkv = nn.Linear(dim, 3 * dim)
        self.out = nn.Linear(dim, dim)
        self.heads, self.p = heads, p

    def forward(self, x):  # x (n, dim); residual pre-norm attention as in TransMIL's TransLayer
        n, dim = x.shape
        q, k, v = self.qkv(self.norm(x)).view(n, 3, self.heads, dim // self.heads).permute(1, 2, 0, 3).unsqueeze(1)
        o = F.scaled_dot_product_attention(q, k, v, dropout_p=self.p if self.training else 0.0)
        return x + self.out(o.squeeze(0).permute(1, 0, 2).reshape(n, dim))


class _PPEG(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.p7 = nn.Conv2d(dim, dim, 7, 1, 3, groups=dim)
        self.p5 = nn.Conv2d(dim, dim, 5, 1, 2, groups=dim)
        self.p3 = nn.Conv2d(dim, dim, 3, 1, 1, groups=dim)

    def forward(self, x, side):
        cls, f = x[:1], x[1:]
        f = f.transpose(0, 1).reshape(1, -1, side, side)
        f = f + self.p7(f) + self.p5(f) + self.p3(f)
        return torch.cat([cls, f.flatten(2).squeeze(0).transpose(0, 1)])


class TransMIL(MIL):
    def __init__(self, d, h=512, heads=8):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU())
        self.cls_token = nn.Parameter(torch.randn(1, h) * 0.02)
        self.l1, self.l2 = _SelfAttention(h, heads), _SelfAttention(h, heads)
        self.ppeg = _PPEG(h)
        self.norm = nn.LayerNorm(h)
        self.head = nn.Linear(h, 1)

    def forward(self, x):
        hdn = self.fc(x)
        n = hdn.shape[0]
        side = math.ceil(math.sqrt(n))
        hdn = torch.cat([hdn, hdn[torch.arange(side * side - n, device=hdn.device) % n]])  # pad to a square
        hdn = torch.cat([self.cls_token.to(hdn.dtype), hdn])
        hdn = self.l2(self.ppeg(self.l1(hdn), side))
        return self.head(self.norm(hdn)[0]).squeeze(), {}


class DSMIL(MIL):
    def __init__(self, d, h=256, q=128, p=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
        self.inst = nn.Linear(h, 1)
        self.q = nn.Sequential(nn.Linear(h, q), nn.ReLU(), nn.Linear(q, q), nn.Tanh())
        self.v = nn.Sequential(nn.Dropout(p), nn.Linear(h, h), nn.ReLU())
        self.bag = nn.Linear(h, 1)
        self.dq = q

    def forward(self, x):
        hdn = self.fc(x)
        c = self.inst(hdn).squeeze(-1)
        Q = self.q(hdn)
        A = torch.softmax(Q @ Q[c.argmax()] / math.sqrt(self.dq), dim=0)
        bag = self.bag((A[:, None] * self.v(hdn)).sum(0)).squeeze()
        mx = c.max()
        return 0.5 * (bag + mx), {"attention": A, "instance": c, "bag": bag, "max": mx}

    def objective(self, logit, aux, y, lossf):
        return 0.5 * lossf(aux["bag"], y) + 0.5 * lossf(aux["max"], y)


class ACMIL(MIL):
    def __init__(self, d, h=256, att=128, p=0.25, branches=5, topk=10, p_mask=0.6, w_div=0.1):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
        self.a = nn.Sequential(nn.Linear(h, att), nn.Tanh())
        self.b = nn.Sequential(nn.Linear(h, att), nn.Sigmoid())
        self.c = nn.Linear(att, branches)
        self.branch = nn.ModuleList([nn.Linear(h, 1) for _ in range(branches)])
        self.cls = nn.Linear(h, 1)
        self.topk, self.p_mask, self.w_div = topk, p_mask, w_div

    def forward(self, x):
        hdn = self.fc(x)
        s = self.c(self.a(hdn) * self.b(hdn)).transpose(0, 1)  # (branches, n)
        if self.training:  # stochastic top-k instance masking
            idx = torch.topk(s, min(self.topk, s.shape[1]), dim=1).indices
            drop = torch.zeros_like(s, dtype=torch.bool)
            drop.scatter_(1, idx, torch.rand(idx.shape, device=s.device) < self.p_mask)
            s = s.masked_fill(drop, -1e4)
        w = torch.softmax(s.float(), dim=1).to(hdn.dtype)
        branch_logits = torch.stack([m(f) for m, f in zip(self.branch, w @ hdn)]).squeeze(-1)
        wm = w.mean(0)
        return self.cls(wm @ hdn).squeeze(), {"attention": wm, "branch_attention": w, "branch_logits": branch_logits}

    def objective(self, logit, aux, y, lossf):
        w = F.normalize(aux["branch_attention"].float(), dim=1)
        sim = w @ w.T
        m = sim.shape[0]
        div = (sim.sum() - sim.diagonal().sum()) / (m * (m - 1))
        branch = torch.stack([lossf(b, y) for b in aux["branch_logits"]]).mean()
        return lossf(logit, y) + branch + self.w_div * div


MODELS = {"abmil": ABMIL, "clam": CLAMSB, "transmil": TransMIL, "dsmil": DSMIL, "acmil": ACMIL}
