"""Do the mutations differ between continents (or countries)?

For every mutation (position + new base) carried by enough genomes: how often it occurs in each group, an overall chi-square test across the groups, and the same test **inside lineages**: the group
labels are shuffled among genomes of the same lineage (a permutation test), so a mutation only counts as different between groups when it differs more than the lineage mix explains. The p values are
corrected for the number of mutations tested (Benjamini-Hochberg). The numbers describe the genomes that were sampled: where sequencing is heavy or early, a group looks different because of who was
sequenced, a mutation frequent in one continent can have been introduced there once and then sampled many times, and a group with few genomes cannot show a difference. They say nothing about where a
mutation arose."""
import math
from collections import Counter, defaultdict

import numpy as np

from . import geo, insight, mutprofile


def _gammaq(a: float, x: float) -> float:
    """Regularised upper incomplete gamma Q(a, x): the chi-square survival function is Q(df / 2, chi2 / 2)."""
    if x <= 0:
        return 1.0
    if x < a + 1:
        s, term, n = 1.0 / a, 1.0 / a, a
        for _ in range(400):
            n += 1
            term *= x / n
            s += term
            if abs(term) < abs(s) * 1e-14:
                break
        return max(0.0, 1.0 - s * math.exp(-x + a * math.log(x) - math.lgamma(a)))
    b, c, d = x + 1 - a, 1e300, 1 / (x + 1 - a)
    h = d
    for i in range(1, 400):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-14:
            break
    return min(1.0, h * math.exp(-x + a * math.log(x) - math.lgamma(a)))


def chi2_sf(chi: float, df: int) -> float:
    return _gammaq(df / 2.0, chi / 2.0)


def bh(ps: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p values (q values), in the order given."""
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    q = [1.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, ps[i] * m / rank)
        q[i] = prev
    return q


def _chi(counts: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    """Chi-square of carriers per group against the group sizes, for the last axis = groups. counts: (..., G), sizes: (G,)."""
    tot = counts.sum(-1, keepdims=True)
    n = sizes.sum()
    exp1 = tot * sizes / n
    exp0 = (n - tot) * sizes / n
    with np.errstate(divide="ignore", invalid="ignore"):
        a = np.where(exp1 > 0, (counts - exp1) ** 2 / exp1, 0.0)
        b = np.where(exp0 > 0, ((sizes - counts) - exp0) ** 2 / exp0, 0.0)
    return (a + b).sum(-1)


def group_of(country: str, by: str) -> str:
    if by == "country":
        pt = geo.resolve(country)
        return pt["name"] if pt else ""
    pt = geo.resolve(country)
    c = pt["continent"] if pt else ""
    return c if c and not c.startswith("Seven seas") else ""


def compare(run: dict, by: str = "continent", min_group: int = 10, min_carriers: int = 5, max_n: float | None = 0.05, adjust: str = "lineage", perms: int = 400, top: int = 80, seed: int = 1) -> dict:
    by = by if by in ("continent", "country") else "continent"
    meta = mutprofile._run_meta(run)
    rows = insight._mutations(run, max_n)
    judged = {s for s in meta}
    ok = mutprofile.good_samples(run, max_n)
    if ok is not None:
        judged &= ok
    grp = {s: group_of(meta[s]["country"], by) for s in judged}
    no_group = sum(1 for s in judged if not grp[s])
    sizes_all = Counter(g for g in grp.values() if g)
    groups = sorted((g for g, n in sizes_all.items() if n >= min_group), key=lambda g: -sizes_all[g])
    small = {g: n for g, n in sizes_all.items() if n < min_group}
    if len(groups) < 2:
        return {"error": f"비교할 수 있는 그룹이 2개 미만입니다(그룹당 유전체 {min_group}건 이상 필요). 그룹을 '나라'로 바꾸거나 최소 유전체 수를 낮춰 보세요.", "sizes": dict(sizes_all)}
    samples = sorted(s for s in judged if grp[s] in groups)
    idx = {s: i for i, s in enumerate(samples)}
    gi = np.array([groups.index(grp[s]) for s in samples])
    G, n = len(groups), len(samples)
    sizes = np.bincount(gi, minlength=G).astype(float)
    sites: dict[tuple, set] = defaultdict(set)
    info: dict[tuple, dict] = {}
    for r in rows:
        i = idx.get(r["sample"])
        if i is None:
            continue
        k = (int(r["pos"]), r["alt"])
        sites[k].add(i)
        info.setdefault(k, r)
    keys = [k for k, v in sites.items() if min_carriers <= len(v) <= n - min_carriers]
    keys.sort()
    if not keys:
        return {"error": "비교할 변이가 없습니다(최소 보유 유전체 수를 낮춰 보세요).", "sizes": dict(sizes_all)}
    X = np.zeros((len(keys), n), dtype=np.float32)
    for j, k in enumerate(keys):
        X[j, list(sites[k])] = 1.0
    onehot = np.eye(G, dtype=np.float32)[gi]                                  # n x G
    counts = X @ onehot                                                       # m x G
    chi = _chi(counts.astype(float), sizes)
    p_raw = [chi2_sf(float(c), G - 1) for c in chi]
    lin = [meta[s]["lineage"] or "(none)" for s in samples]
    p_adj = None
    if adjust == "lineage":
        rng = np.random.default_rng(seed)
        strata: dict[str, list] = defaultdict(list)
        for i, l in enumerate(lin):
            strata[l].append(i)
        ge = np.zeros(len(keys))
        done = 0
        for _ in range(perms):
            lab = gi.copy()
            for ix in strata.values():
                if len(ix) > 1:
                    lab[ix] = gi[rng.permutation(ix)]
            c = X @ np.eye(G, dtype=np.float32)[lab]
            ge += _chi(c.astype(float), sizes) >= chi - 1e-9
            done += 1
        p_adj = (ge + 1) / (done + 1)
    q_raw = bh(p_raw)
    q_adj = bh([float(x) for x in p_adj]) if p_adj is not None else None
    out = []
    for j, k in enumerate(keys):
        r = info[k]
        c = counts[j]
        fr = c / sizes
        out.append({"key": f"{k[0]}{k[1]}", "label": f"{r['ref']}{k[0]}{k[1]}", "pos": k[0], "ref": r["ref"], "alt": k[1], "gene": r["gene"], "aa": r["aa_change"], "effect": r["effect"],
                    "apobec": int(r["profile_of_interest"] not in ("", "other")), "n": int(c.sum()), "counts": [int(x) for x in c], "freq": [round(float(x), 4) for x in fr],
                    "chi2": round(float(chi[j]), 2), "p": p_raw[j], "q": q_raw[j], "p_adj": float(p_adj[j]) if p_adj is not None else None, "q_adj": q_adj[j] if q_adj is not None else None,
                    "range": round(float(fr.max() - fr.min()), 4), "top": groups[int(fr.argmax())]})
    key = "q_adj" if adjust == "lineage" else "q"
    out.sort(key=lambda x: (x[key], -x["range"]))
    nsig = sum(1 for x in out if x["q"] < 0.05)
    nsig_adj = sum(1 for x in out if x["q_adj"] is not None and x["q_adj"] < 0.05)
    per = []
    for g in groups:
        ix = [i for i, s in enumerate(samples) if grp[s] == g]
        nm = np.asarray(X[:, ix].sum(0)).ravel() if ix else np.zeros(0)
        lc = Counter(lin[i] for i in ix).most_common(4)
        per.append({"group": g, "n": len(ix), "mut_median": float(np.median(nm)) if len(nm) else 0.0, "lineages": [[a, b] for a, b in lc],
                    "countries": [c for c, _ in Counter(meta[samples[i]]["country"] for i in ix).most_common(3)]})
    return {"by": by, "groups": groups, "sizes": [int(x) for x in sizes], "per_group": per, "genomes": n, "no_group": no_group, "small_groups": small, "tested": len(out), "adjust": adjust,
            "perms": perms if adjust == "lineage" else 0, "significant_raw": nsig, "significant_adjusted": nsig_adj if adjust == "lineage" else None,
            "rows": out[:top], "n_lineages": len(set(lin)), "min_carriers": min_carriers}
