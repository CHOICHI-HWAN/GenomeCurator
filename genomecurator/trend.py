"""Is a mutation rising among the sampled genomes, once the country mix is held fixed? (evidence E5 of docs/ADAPTIVE.md)

Logistic regression of "carries the mutation" on time (months) with one intercept per country (a mild ridge keeps a country where nobody / everybody carries it finite). The slope is the change of
the log-odds per month with a Wald 95% interval; an interval above zero means the share rises more than the changing mix of countries explains. It cannot correct for sampling that changes inside a
country (a hospital that starts sequencing), and the model assumes a straight line in the log-odds."""
import math

import numpy as np


def context(samples: list[str], month_of: dict[str, int], country_of: dict[str, str]) -> dict:
    """Design matrix pieces shared by all mutations of one pool: month index (centred) and country dummies."""
    ids = [s for s in samples if s in month_of]
    t = np.array([month_of[s] for s in ids], dtype=float)
    t -= t.mean() if len(t) else 0.0
    cs = sorted({country_of.get(s) or "(없음)" for s in ids})
    D = np.zeros((len(ids), len(cs)))
    for i, s in enumerate(ids):
        D[i, cs.index(country_of.get(s) or "(없음)")] = 1.0
    return {"ids": ids, "index": {s: i for i, s in enumerate(ids)}, "t": t, "D": D}


def carrier_trend(carriers: set[str], ctx: dict, ridge: float = 1e-2) -> dict | None:
    """{slope, lo, hi, n} (log-odds per month) or None when there are too few carriers or only one month."""
    n = len(ctx["ids"])
    y = np.zeros(n)
    for s in carriers:
        i = ctx["index"].get(s)
        if i is not None:
            y[i] = 1.0
    if y.sum() < 3 or y.sum() > n - 3 or np.ptp(ctx["t"]) == 0:
        return None
    X = np.concatenate([ctx["t"][:, None], ctx["D"]], axis=1)
    pen = np.full(X.shape[1], ridge)
    pen[0] = 0.0
    beta = np.zeros(X.shape[1])
    for _ in range(60):
        eta = np.clip(X @ beta, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        W = p * (1 - p) + 1e-9
        g = X.T @ (y - p) - pen * beta
        H = (X.T * W) @ X + np.diag(pen)
        step = np.linalg.solve(H, g)
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break
    cov = np.linalg.inv(H)
    se = math.sqrt(max(0.0, cov[0, 0]))
    return {"slope": float(beta[0]), "lo": float(beta[0] - 1.96 * se), "hi": float(beta[0] + 1.96 * se), "n": n}
