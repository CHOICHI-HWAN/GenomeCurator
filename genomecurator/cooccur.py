"""Which mutations appear together: groups of co-occurring mutations and a principal-component view of the genomes.

The input is the run's genome x mutation table (1 = the genome carries the mutation; ITRs left out, genomes with at most `max_n` N). Two views of the same matrix:

* **Mutation groups** - mutations are clustered by how much their carrier sets overlap (Jaccard distance, average linkage, cut at `merge`): mutations on the same branch of the tree
  have the same carriers and fall into one group, mutations that mostly travel together (a few genomes lack one of them, e.g. for N) join it. This is the exact answer to "which
  mutations come together"; PCA would blur the partial overlaps.
* **Principal components of the genomes** - PCA (singular value decomposition of the centred 0/1 table) places each genome by the mutations it carries; clades show up as clumps, the
  loadings say which mutations make an axis, and a mutation that repeats independently of the tree (typically APOBEC3 hot spots) stands off the lineage axes.

What it does not say: co-occurring means *on the same branch*, not linked in function or under selection. Rare mutations (fewer than `min_carriers` genomes) are left out - they are noise
for both views. A genome with many N's has unknown, not absent, mutations, hence the N limit."""
from collections import Counter, defaultdict

import numpy as np

from . import insight, lineages, log, mutprofile, novel, outbreak

_log = log.get("cooccur")


def _upgma(dist: np.ndarray, cut: float) -> list[list[int]]:
    """Average-linkage clusters of the items, merging while the closest pair is at most `cut` apart."""
    m = dist.shape[0]
    clusters = [[i] for i in range(m)]
    d = dist.astype(np.float64).copy()
    np.fill_diagonal(d, np.inf)
    active = list(range(m))
    size = {i: 1 for i in range(m)}
    while len(active) > 1:
        sub = d[np.ix_(active, active)]
        k = int(np.argmin(sub))
        a, b = divmod(k, len(active))
        if sub[a, b] > cut:
            break
        ia, ib = active[a], active[b]
        for j in active:
            if j not in (ia, ib):
                d[ia, j] = d[j, ia] = (d[ia, j] * size[ia] + d[ib, j] * size[ib]) / (size[ia] + size[ib])
        size[ia] += size[ib]
        clusters[ia] += clusters[ib]
        clusters[ib] = []
        active.remove(ib)
        d[ib, :] = d[:, ib] = np.inf
    return [clusters[i] for i in active]


def analyse(run: dict, max_n: float | None = 0.05, min_carriers: int = 3, max_share: float = 0.9, merge: float = 0.15, max_sites: int = 400) -> dict:
    rows = insight._mutations(run, max_n)
    samples = sorted({r["sample"] for r in rows})
    n = len(samples)
    if n < 5:
        return {"available": False, "reason": "분석할 유전체가 5건 미만입니다 (N 비율 상한을 풀어 보세요)."}
    si = {s: i for i, s in enumerate(samples)}
    carriers: dict[tuple, set] = defaultdict(set)
    info: dict[tuple, dict] = {}
    for r in rows:
        k = (int(r["pos"]), r["alt"])
        carriers[k].add(r["sample"])
        info.setdefault(k, r)
    sites = [k for k, v in carriers.items() if min_carriers <= len(v) <= max_share * n]
    sites.sort(key=lambda k: -len(carriers[k]))
    sites = sites[:max_sites]
    if len(sites) < 3:
        return {"available": False, "reason": "묶을 변이가 3개 미만입니다 (최소 보유 수를 낮추거나 N 비율 상한을 풀어 보세요)."}
    m = len(sites)
    X = np.zeros((n, m), dtype=np.float32)
    for j, k in enumerate(sites):
        for s in carriers[k]:
            X[si[s], j] = 1
    inter = X.T @ X
    sz = np.diag(inter)
    union = sz[:, None] + sz[None, :] - inter
    jac = np.where(union > 0, inter / np.maximum(union, 1), 0.0)
    groups_idx = _upgma(1.0 - jac, merge)
    # PCA of the genomes
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    ev = (S ** 2) / max(np.sum(S ** 2), 1e-12)
    kpc = min(6, len(S))
    scores = U[:, :kpc] * S[:kpc]
    load = Vt[:kpc].T                                                           # sites x components
    gs = {g["id"]: g for g in lineages._genomes(run)}
    ref = mutprofile.summary(run).get("reference", "")
    cand_by_label: dict[str, list[str]] = {}
    try:                                                                         # link to the new-variant candidates (labels are comparable when the tree uses the run's reference)
        nd = novel.discover(run, 120, 3, max_n)
        if nd.get("available") and nd.get("reference") == ref:
            for c in nd["candidates"]:
                for mm in c["mutations"]:
                    cand_by_label.setdefault(mm["label"], []).append(c["id"])
    except Exception:
        _log.warning("new-variant candidates unavailable for the PCA page", exc_info=True)
    label = lambda k: f"{info[k]['ref']}{k[0]}{k[1]}"
    out_groups = []
    order = []
    for gi, idxs in enumerate(sorted(groups_idx, key=lambda ix: (-len(carriers[sites[ix[0]]]), sites[ix[0]][0]))):
        ks = [sites[i] for i in sorted(idxs, key=lambda i: sites[i][0])]
        cs = set.union(*(carriers[k] for k in ks))
        core = set.intersection(*(carriers[k] for k in ks))
        meta = [gs[s] for s in cs if s in gs]
        ds = sorted(d for d in (lineages.parse_date(g["date"])[3] for g in meta) if d)
        cl = [cand for k in ks for cand in cand_by_label.get(label(k), [])]
        out_groups.append({"id": f"G{gi + 1}", "mutations": [{"label": label(k), "gene": info[k]["gene"], "aa": info[k]["aa_change"], "effect": info[k]["effect"], "apobec": insight._is_apo(info[k]), "n": len(carriers[k])} for k in ks],
                           "n_mutations": len(ks), "n_genomes": len(cs), "n_all_carry": len(core), "similarity": round(float(np.mean([jac[sites.index(a), sites.index(b)] for a in ks for b in ks if a != b])), 3) if len(ks) > 1 else 1.0,
                           "lineages": dict(Counter(g["lineage"] for g in meta if g["lineage"]).most_common(4)), "countries": dict(Counter(g["country"] for g in meta if g["country"]).most_common(5)),
                           "first": ds[0] if ds else "", "last": ds[-1] if ds else "", "apobec": sum(1 for k in ks if insight._is_apo(info[k])), "candidates": sorted(set(cl))})
        order += [sites.index(k) for k in ks]
    gid = {mm["label"]: g["id"] for g in out_groups for mm in g["mutations"]}
    muts = [{"label": label(k), "gene": info[k]["gene"], "aa": info[k]["aa_change"], "apobec": insight._is_apo(info[k]), "n": len(carriers[k]), "group": gid[label(k)],
             "loading": [round(float(load[j, c]), 4) for c in range(min(4, kpc))]} for j, k in enumerate(sites)]
    lin_of = [gs.get(s, {}).get("lineage", "") for s in samples]
    genomes = [{"id": s, "lineage": lin_of[i], "country": gs.get(s, {}).get("country", ""), "date": gs.get(s, {}).get("date", ""), "pc": [round(float(scores[i, c]), 3) for c in range(min(4, kpc))],
                "n_mut": int(X[i].sum())} for i, s in enumerate(samples)]
    # variance explained by lineage per component
    eta = []
    for c in range(min(4, kpc)):
        pc = scores[:, c]
        tot = float(((pc - pc.mean()) ** 2).sum()) or 1.0
        within = 0.0
        for L in set(lin_of):
            v = pc[[i for i, l in enumerate(lin_of) if l == L]]
            within += float(((v - v.mean()) ** 2).sum())
        eta.append(round(1 - within / tot, 3))
    heat = {"order": [label(sites[j]) for j in order], "groups": [gid[label(sites[j])] for j in order],
            "values": [[round(float(jac[a, b]), 2) for b in order] for a in order] if m <= 200 else None}
    return {"available": True, "reference": ref, "n_genomes": n, "n_sites": m, "max_n": max_n, "min_carriers": min_carriers, "merge": merge, "explained": [round(float(x), 4) for x in ev[:kpc]],
            "lineage_eta2": eta, "groups": out_groups, "mutations": muts, "genomes": genomes, "heat": heat,
            "identical_groups": sum(1 for g in out_groups if g["n_mutations"] > 1), "grouped_mutations": sum(g["n_mutations"] for g in out_groups if g["n_mutations"] > 1)}
