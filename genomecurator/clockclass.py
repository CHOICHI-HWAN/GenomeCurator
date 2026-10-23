"""Does the molecular clock tick in APOBEC3-type changes, in the other changes, or in both?

In human mpox transmission most changes are made by the host's APOBEC3 enzymes (TC>TT, GA>AA), at a rate that is roughly constant, and the rest are rare. A clock fitted to all changes together mixes the two:
the APOBEC3-type changes carry the temporal signal, the others add noise (or nothing). This module counts, for every dated tip of the dataset's Nextstrain tree, the APOBEC3-type and the other substitutions on its
path from the root, and regresses each count on the collection date - the root-to-tip regression of TempEst, once for each class (and for the lineages with enough tips).

What it shows: the slope (changes per year per genome), R^2, the date where the line meets zero (an estimate of when the counted changes began to accumulate - the start of APOBEC3-driven human
transmission when the APOBEC3 line is the clean one), and how many changes there are at all. What it does not: a regression through tips is not a tree-aware estimate (tips of one lineage share their path, so the
intervals are too narrow); the date where a line meets zero is not a tMRCA; APOBEC3-type mutations recur at favoured sites, so they saturate over longer times; and the start of the counting depends on the tree's root."""
import sys
from collections import Counter
from datetime import date

import numpy as np

from . import insight, lineages, outbreak


def _tip_rows(run: dict) -> list[dict]:
    tree, root, _ns = outbreak._tree(run)
    if not tree:
        return []
    meta = {g["id"]: g for g in lineages._genomes(run)}
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    rows = []

    def walk(n, apo, oth):
        for ref, pos, alt in outbreak._nuc_muts(n):
            if insight._apo_branch(ref, alt, pos, root):
                apo += 1
            else:
                oth += 1
        kids = n.get("children") or []
        if not kids:
            g = meta.get(n["name"])
            if g:
                kind, _wk, _mo, shown = lineages.parse_date(g["date"])
                if kind == "day":
                    y, m, d = (int(x) for x in shown.split("-"))
                    rows.append({"id": n["name"], "t": y + (date(y, m, d).timetuple().tm_yday - 0.5) / 365.25, "apobec": apo, "other": oth, "lineage": g["lineage"], "country": g["country"]})
        for c in kids:
            walk(c, apo, oth)
    walk(tree, 0, 0)
    return rows


def regress(x: np.ndarray, y: np.ndarray) -> dict:
    n = len(x)
    if n < 5 or np.ptp(x) == 0:
        return {"n": n, "slope": None, "r2": None, "lo": None, "hi": None, "zero_at": None, "mean": float(y.mean()) if n else None}
    A = np.vstack([x, np.ones(n)]).T
    b, c = np.linalg.lstsq(A, y, rcond=None)[0]
    res = y - (b * x + c)
    se = float(np.sqrt(res.var(ddof=2) / ((x - x.mean()) ** 2).sum())) if n > 2 else 0.0
    r2 = float(np.corrcoef(x, y)[0, 1] ** 2) if y.std() > 0 else 0.0
    return {"n": n, "slope": float(b), "r2": r2, "lo": float(b - 1.96 * se), "hi": float(b + 1.96 * se), "zero_at": float(-c / b) if abs(b) > 1e-9 else None, "mean": float(y.mean())}


def analyze(run: dict, min_lineage_tips: int = 30) -> dict:
    rows = _tip_rows(run)
    if not rows:
        return {"available": False, "reason": "이 데이터셋의 Nextstrain 분석(트리)이 필요합니다."}
    t = np.array([r["t"] for r in rows])
    apo = np.array([r["apobec"] for r in rows], dtype=float)
    oth = np.array([r["other"] for r in rows], dtype=float)
    out = {"available": True, "tips": len(rows), "first": float(t.min()), "last": float(t.max()), "classes": {"apobec": regress(t, apo), "other": regress(t, oth), "all": regress(t, apo + oth)},
           "apobec_share": float(apo.sum() / max(1.0, apo.sum() + oth.sum()))}
    by = Counter(r["lineage"] for r in rows)
    out["lineages"] = []
    for lin, n in by.most_common():
        if n >= min_lineage_tips and lin:
            ix = [i for i, r in enumerate(rows) if r["lineage"] == lin]
            out["lineages"].append({"lineage": lin, "apobec": regress(t[ix], apo[ix]), "other": regress(t[ix], oth[ix])})
    step = max(1, len(rows) // 400)
    out["points"] = [[round(r["t"], 3), int(r["apobec"]), int(r["other"])] for r in rows[::step]]
    return out
