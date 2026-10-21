"""The clade / lineage hierarchy of mpox as a diagram, and the growth of lineages against each other.

* `tree(run)` - the hierarchy mpox > clade (I: Ia, Ib; II: IIa, IIb) > outbreak (the lineage-designation outbreaks sh2017 / sh2023 / sh2024) > lineage (parent links from the mpxv-lineages definitions),
  with, for every node, the mutations that *introduced* it (new relative to its parent), the cumulative number of defining mutations along the path from the root, and the genomes of the run (own,
  with all sub-lineages, and within a recent window). That is what a tree diagram, a Sankey diagram (link width = genomes) and the accumulation of mutations along a path are drawn from.
* `growth(run, ...)` - a multinomial logistic regression of the lineage shares over time (the method of the UK and Nextstrain growth-advantage estimates): for each lineage the log-odds against a reference
  lineage change linearly with time. Reported: the slope against the reference and against the average lineage (per week or month) with 95% intervals, and the fitted and observed shares.

What the growth numbers are *not*: they describe how the shares in the **sample** change. Sampling is uneven in time and place and not random, a share rises when sequencing of one place rises, the model
assumes a constant advantage over the window, and a lineage with few genomes has a wide interval. They are not transmissibility estimates."""
from collections import Counter, defaultdict

import numpy as np

from . import lineages, mutprofile

CLADES = {"I": ["Ia", "Ib"], "II": ["IIa", "IIb"]}
_CLADE_PARENT = {c: g for g, cs in CLADES.items() for c in cs}


def _node(id_, label, kind, parent, clade="", outbreak="", date=""):
    return {"id": id_, "label": label, "kind": kind, "parent": parent, "clade": clade, "outbreak": outbreak, "designation": date, "own_snps": [], "cum_snps": 0, "n_own": 0, "n_sub": 0, "n_win": 0,
            "n_sub_win": 0, "first": "", "last": "", "defined": True}


def _skeleton() -> dict[str, dict]:
    nodes = {"mpox": _node("mpox", "Mpox virus", "root", "")}
    for g, cs in CLADES.items():
        nodes[g] = _node(g, f"Clade {g}", "clade", "mpox", clade=g)
        for c in cs:
            nodes[c] = _node(c, f"Clade {c}", "clade", g, clade=c)
    return nodes


def tree(run: dict, weeks: int = 0) -> dict:
    """The hierarchy with genome counts. `weeks` > 0: also count the genomes collected in the last that many weeks before the newest collection date of the run (0: no window)."""
    defs = lineages.definitions()
    nodes = _skeleton()
    try:
        _name, ncseq, cds = lineages.nc_reference()
    except RuntimeError:
        ncseq, cds = "", []
    for ob, d in sorted(defs.items()):
        root_def = d.get("", {})
        clade = root_def.get("clade") or next((x["clade"] for x in d.values() if x["clade"]), "")
        rid = f"ob:{ob}"
        nodes[rid] = _node(rid, ob, "outbreak", clade if clade in nodes else "mpox", clade=clade, outbreak=ob, date=root_def.get("designation_date", ""))
        for name, x in d.items():
            if name:
                nodes[f"{ob}/{name}"] = _node(f"{ob}/{name}", name, "lineage", "", clade=x["clade"], outbreak=ob, date=x["designation_date"])
        for name, x in d.items():
            if not name:
                continue
            par = x["parent"] if x["parent"] not in ("None", "none", "null") else ""
            nodes[f"{ob}/{name}"]["parent"] = f"{ob}/{par}" if par else rid
            if par and par not in d and f"{ob}/{par}" not in nodes:            # a parent without a definition of its own (an alias root such as B): kept as a node so the hierarchy is connected
                nodes[f"{ob}/{par}"] = _node(f"{ob}/{par}", par, "lineage", rid, clade=x["clade"], outbreak=ob)
                nodes[f"{ob}/{par}"]["defined"] = False
        for name, x in d.items():
            nid = f"{ob}/{name}" if name else rid
            snps = lineages.defining_snps(ob, name) if name else [(p, n, ob) for p, n in x["snps"]]
            nodes[nid]["cum_snps"] = len(snps)
            for pos, nuc, by in snps:
                if (by or ob) == (name or ob):
                    ref_b = ncseq[pos - 1] if 0 < pos <= len(ncseq) else "?"
                    an = mutprofile.annotate(pos, ref_b, nuc, cds, ncseq) if ref_b in "ACGT" else {"gene": "", "aa_change": "", "effect": ""}
                    nodes[nid]["own_snps"].append({"label": f"{ref_b}{pos}{nuc}", "gene": an["gene"], "aa": an["aa_change"], "effect": an["effect"]})
    for n in nodes.values():                                            # implicit nodes get the cumulative count of their parent chain
        if n["kind"] == "lineage" and not n["defined"]:
            n["cum_snps"] = 0
    gs = lineages._genomes(run)
    dated = [(g, lineages.parse_date(g["date"])) for g in gs]
    days = sorted(p[3] for _g, p in dated if p[0] == "day")
    cutoff = ""
    if weeks and days:
        from datetime import date, timedelta
        y, m, d = (int(x) for x in days[-1].split("-"))
        cutoff = (date(y, m, d) - timedelta(days=7 * weeks)).isoformat()
    for g, p in dated:
        nid = f"{g['outbreak']}/{g['lineage']}" if g["lineage"] and f"{g['outbreak']}/{g['lineage']}" in nodes else (f"ob:{g['outbreak']}" if f"ob:{g['outbreak']}" in nodes else (g["clade"] if g["clade"] in nodes else "mpox"))
        n = nodes[nid]
        n["n_own"] += 1
        if p[0]:
            n["first"] = min(filter(None, [n["first"], p[3]]))
            n["last"] = max(n["last"], p[3])
            if cutoff and p[0] == "day" and p[3] >= cutoff:
                n["n_win"] += 1
    kids: dict[str, list] = defaultdict(list)
    for n in nodes.values():
        if n["parent"]:
            kids[n["parent"]].append(n["id"])

    def sub(i, seen=()):
        n = nodes[i]
        n["n_sub"] = n["n_own"] + sum(sub(k, seen + (i,)) for k in kids[i] if k not in seen)
        n["n_sub_win"] = n["n_win"] + sum(nodes[k]["n_sub_win"] for k in kids[i] if k not in seen)
        return n["n_sub"]
    sub("mpox")
    order = []

    def walk(i, depth):
        nodes[i]["depth"] = depth
        order.append(i)
        for k in sorted(kids[i], key=lambda k: (-nodes[k]["n_sub"], nodes[k]["label"])):
            walk(k, depth + 1)
    walk("mpox", 0)
    for i in nodes:
        if i not in order:                                              # defensive: a cycle or an orphan in the definitions
            nodes[i]["depth"], nodes[i]["parent"] = 1, "mpox"
            order.append(i)
    keep = [n for n in (nodes[i] for i in order) if n["kind"] in ("root", "clade") or n["n_sub"] or n["kind"] != "lineage" or True]
    return {"nodes": keep, "genomes": len(gs), "dated": sum(1 for _g, p in dated if p[0]), "window_weeks": weeks, "cutoff": cutoff, "definitions": bool(defs), "commit": lineages.commit() if defs else ""}


# ------------------------------------------------------------------------------------------------ growth advantage
def _fit(counts: np.ndarray, t: np.ndarray, ridge: float = 1e-4):
    """Multinomial logistic regression, log(p_k / p_0) = a_k + r_k * t. counts: periods x classes (class 0 = reference). -> (theta (K-1, 2), covariance of theta.ravel())."""
    T, K = counts.shape
    X = np.stack([np.ones(T), t], axis=1)
    theta = np.zeros((K - 1, 2))
    N = counts.sum(1)
    for _ in range(60):
        eta = np.concatenate([np.zeros((T, 1)), X @ theta.T], axis=1)
        eta -= eta.max(1, keepdims=True)
        p = np.exp(eta)
        p /= p.sum(1, keepdims=True)
        grad = ((counts[:, 1:] - N[:, None] * p[:, 1:]).T @ X).ravel() - ridge * theta.ravel()
        H = np.zeros(((K - 1) * 2, (K - 1) * 2))
        for k in range(1, K):
            for l in range(1, K):
                w = N * ((p[:, k] if k == l else 0) - p[:, k] * p[:, l])
                H[(k - 1) * 2:k * 2, (l - 1) * 2:l * 2] = -(X.T * w) @ X
        H -= ridge * np.eye(len(H))
        step = np.linalg.solve(H, -grad)
        theta += step.reshape(theta.shape)
        if np.max(np.abs(step)) < 1e-8:
            break
    cov = np.linalg.inv(-H)
    return theta, cov


def growth(run: dict, bin_: str = "week", weeks: int = 26, country: str = "", min_n: int = 20, min_periods: int = 3, ref: str = "") -> dict:
    """Growth advantage of the lineages of a run over the last `weeks` weeks (months when bin_ = 'month': the same span in months). Lineages with fewer than `min_n` genomes in the window or fewer than
    `min_periods` occupied periods are folded into their closest modelled ancestor (else into 'other'), as the UK estimates do. `country` fits the genomes of one country only."""
    bin_ = bin_ if bin_ in ("week", "month") else "week"
    gs = [g for g in lineages._genomes(run) if (not country or g["country"] == country)]
    rows = []
    for g in gs:
        kind, wk, mo, shown = lineages.parse_date(g["date"])
        key = wk if bin_ == "week" else mo
        if key and g["lineage"] and (bin_ == "month" or kind == "day"):
            rows.append((key, g["lineage"], g["outbreak"], shown))
    if not rows:
        return {"error": "채집일과 lineage가 있는 유전체가 없습니다."}
    last = max(r[0] for r in rows)
    allp = lineages._weeks(min(r[0] for r in rows), last) if bin_ == "week" else lineages._months(min(r[0] for r in rows), last)
    span = weeks if bin_ == "week" else max(2, round(weeks / 4.345))
    periods = allp[-span:]
    pset = set(periods)
    win = [r for r in rows if r[0] in pset]
    tot = Counter(r[1] for r in win)
    occ = defaultdict(set)
    for r in win:
        occ[r[1]].add(r[0])
    defs = lineages.definitions()
    parent = {}
    for ob, d in defs.items():
        for name, x in d.items():
            if name and x["parent"]:
                parent.setdefault(name, x["parent"])
    ok = {l for l, n in tot.items() if n >= min_n and len(occ[l]) >= min_periods}

    def home(l):
        seen = set()
        while l and l not in ok and l not in seen:
            seen.add(l)
            l = parent.get(l, "")
        return l if l in ok else "other"
    grouped = Counter()
    cell: dict[tuple, int] = Counter()
    for key, lin, _ob, _s in win:
        h = home(lin)
        cell[(h, key)] += 1
        grouped[h] += 1
    classes = [c for c in sorted(grouped, key=lambda c: -grouped[c]) if c != "other"]
    if len(classes) < 2:
        return {"error": f"이 구간에서 모델에 넣을 수 있는 lineage가 2개 미만입니다(각각 {min_n}건 이상, {min_periods}개 기간 이상 필요). 기간을 늘리거나 최소 건수를 낮춰 보세요.", "periods": periods, "classes": classes}
    if grouped.get("other"):
        classes.append("other")
    refc = ref if ref in classes else classes[0]
    classes = [refc] + [c for c in classes if c != refc]
    used = [p for p in periods if any(cell[(c, p)] for c in classes)]
    C = np.array([[cell[(c, p)] for c in classes] for p in used], dtype=float)
    t = np.arange(len(used), dtype=float) - (len(used) - 1)                # 0 at the last period
    theta, cov = _fit(C, t)
    K = len(classes)
    eta = np.concatenate([np.zeros((len(used), 1)), np.stack([np.ones(len(used)), t], 1) @ theta.T], 1)
    eta -= eta.max(1, keepdims=True)
    P = np.exp(eta)
    P /= P.sum(1, keepdims=True)
    pl = P[-1]
    r = np.concatenate([[0.0], theta[:, 1]])
    mean_r = float(pl @ r)
    out = []
    for k, c in enumerate(classes):
        if k == 0:
            se_ref = 0.0
        else:
            se_ref = float(np.sqrt(cov[(k - 1) * 2 + 1, (k - 1) * 2 + 1]))
        g = np.zeros((K - 1) * 2)                                           # advantage over the average lineage: r_k - sum_j p_j r_j
        for j in range(1, K):
            g[(j - 1) * 2 + 1] = (1.0 if j == k else 0.0) - pl[j]
        se_mean = float(np.sqrt(max(0.0, g @ cov @ g)))
        d_ref, d_mean = float(r[k]), float(r[k] - mean_r)
        out.append({"lineage": c, "n": int(grouped[c]), "share_last": float(pl[k]), "obs_last": float(C[-1, k] / max(1, C[-1].sum())),
                    "vs_ref": d_ref, "vs_ref_lo": d_ref - 1.96 * se_ref, "vs_ref_hi": d_ref + 1.96 * se_ref,
                    "vs_mean": d_mean, "vs_mean_lo": d_mean - 1.96 * se_mean, "vs_mean_hi": d_mean + 1.96 * se_mean,
                    "pct_per_period": float(np.expm1(d_mean) * 100),
                    "modelled": c != "other", "folded": sorted(l for l in tot if l not in ok and home(l) == c)[:12] if c != "other" else sorted(l for l in tot if l not in ok and home(l) == "other")[:12],
                    "observed": [int(x) for x in C[:, k]], "fitted": [float(x) for x in P[:, k]]})
    warnings = []
    if len(used) < 6:
        warnings.append(f"기간이 {len(used)}개뿐이라 추세를 믿기 어렵습니다.")
    if int(C.sum()) < 200:
        warnings.append(f"모델에 들어간 유전체가 {int(C.sum())}건뿐이라 구간이 넓습니다.")
    return {"bin": bin_, "periods": used, "classes": classes, "reference": refc, "rows": out, "genomes": int(C.sum()), "window": span, "country": country, "min_n": min_n, "warnings": warnings,
            "totals": [int(x) for x in C.sum(1)]}
