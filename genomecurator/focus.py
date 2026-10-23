"""One group of genomes against the rest of the dataset: a country's genomes (or only the private ones) seen from their own side.

For a group (default: the genomes of one country, optionally only the private ones) of a mutation-profile run:
* **profile**: how many mutations each genome has, which of them are APOBEC3-type, and which mutations *all* genomes of the group share;
* **distinct mutations**: every mutation carried by the group, with how often the others carry it - overall and among the others of the **same lineage** (what a lineage's own mutations do not explain) -
  and a class: carried by no other genome / by under 2% of them / by under 2% of the others of the same lineage / common. Mutations that nearly all the others of the group's lineages carry but the
  group lacks are listed as well (a lost marker or missing data);
* **nearest relatives**: for each group genome the genomes with the fewest differences in mutations (symmetric difference of the mutation sets), with country, lineage and date, and how many genomes are *equally* near (a tie: the order among them, by date, says nothing);
* for the few most distinct mutations, whether the callers (the run's alignment, snippy, lofreq, Nextclade) agree, and how often the tree says the mutation arose independently.

Limits: a group of three genomes cannot support a statistical test (no p value is given for groups under ten); "nearest" depends on who was sequenced and is not a transmission route; a mutation seen
in no other genome may simply not have been sampled elsewhere."""
from collections import Counter, defaultdict

from . import geo, insight, lineages, mutprofile, origins, variantcall

MIN_FOR_TEST = 10


def group_members(meta: dict, country: str = "South Korea", private_only: bool = False) -> list[str]:
    want = (geo.resolve(country) or {}).get("name") or country
    out = []
    for s, m in meta.items():
        pt = geo.resolve(m["country"])
        if (pt["name"] if pt else m["country"]) != want:
            continue
        if private_only and not str(s).startswith("private:") and not _is_private(s):
            continue
        out.append(s)
    return sorted(out)


def _is_private(sample: str) -> bool:
    """Public accessions of Pathoplexus / INSDC look like PP_xxxx / PV123456 / PZ123456 ...; the private store's ids are free text (e.g. a lab number)."""
    import re
    return not re.fullmatch(r"(PP_[0-9A-Z]+|[A-Z]{1,2}_?[0-9]{5,8})(\.[0-9]+)?", str(sample))


def analyze(run: dict, country: str = "South Korea", private_only: bool = False, max_n: float | None = None, top: int = 40) -> dict:
    meta = mutprofile._run_meta(run)
    members = group_members(meta, country, private_only)
    if not members:
        return {"error": f"이 실행에 '{country}'{' 비공개' if private_only else ''} 유전체가 없습니다.", "countries": _countries(meta)}
    rows = insight._mutations(run, max_n)
    ok = mutprofile.good_samples(run, max_n)
    pool = [s for s in meta if ok is None or s in ok or s in members]
    G = set(members)
    O = [s for s in pool if s not in G]
    sites: dict[tuple, set] = defaultdict(set)
    info: dict[tuple, dict] = {}
    gm: dict[str, set] = defaultdict(set)
    for r in rows:
        k = (int(r["pos"]), r["alt"])
        sites[k].add(r["sample"])
        info.setdefault(k, r)
        gm[r["sample"]].add(k)
    glin = Counter(meta[s]["lineage"] or "(없음)" for s in members)
    same = [s for s in O if (meta[s]["lineage"] or "(없음)") in glin]
    shared = set.intersection(*(gm[s] for s in members)) if members else set()
    apo = lambda k: info[k]["profile_of_interest"] not in ("", "other")
    per = [{"id": s, "lineage": meta[s]["lineage"], "clade": meta[s]["clade"], "date": meta[s]["date"], "mutations": len(gm[s]), "apobec": sum(1 for k in gm[s] if apo(k))} for s in members]
    d_rows = []
    for k in set().union(*(gm[s] for s in members)):
        v = sites[k]
        kc, oc, sc = len(v & G), len(v & set(O)), len(v & set(same))
        x = info[k]
        cls = ("only" if oc == 0 else "rare" if oc / max(1, len(O)) < 0.02 else "rare_in_lineage" if same and sc / len(same) < 0.02 else "common")
        d_rows.append({"key": f"{k[0]}{k[1]}", "label": f"{x['ref']}{k[0]}{k[1]}", "pos": k[0], "ref": x["ref"], "alt": k[1], "gene": x["gene"], "aa": x["aa_change"], "effect": x["effect"], "apobec": int(apo(k)),
                       "in_group": kc, "group_n": len(members), "others": oc, "others_n": len(O), "same_lineage": sc, "same_lineage_n": len(same), "class": cls})
    order = {"only": 0, "rare": 1, "rare_in_lineage": 2, "common": 3}
    d_rows.sort(key=lambda r: (order[r["class"]], r["others"], r["pos"]))
    lost = []
    if same:
        for k, v in sites.items():
            sc = len(v & set(same))
            if sc >= 0.9 * len(same):
                miss = [s for s in members if s not in v]
                if miss:
                    x = info[k]
                    lost.append({"label": f"{x['ref']}{k[0]}{k[1]}", "gene": x["gene"], "aa": x["aa_change"], "missing_in": miss, "same_lineage": sc, "same_lineage_n": len(same)})
    near, ties = {}, {}
    for s in members:
        ds = sorted(((len(gm[s] ^ gm[o]), o) for o in O), key=lambda t: (t[0], str(meta[t[1]]["date"])))[:8]
        near[s] = [{"id": o, "distance": d, "country": meta[o]["country"], "lineage": meta[o]["lineage"], "date": meta[o]["date"]} for d, o in ds]
        alld = [(len(gm[s] ^ gm[o]), o) for o in O]
        best = min((d for d, _o in alld), default=None)
        tied = [o for d, o in alld if d == best]
        ties[s] = {"distance": best, "n": len(tied), "countries": dict(Counter((geo.resolve(meta[o]["country"]) or {}).get("name") or meta[o]["country"] for o in tied).most_common(8)),
                   "first": min((str(meta[o]["date"]) for o in tied), default="")}
    tree = origins.independent_origins(run)
    for r in d_rows[:12]:
        if r["class"] in ("only", "rare", "rare_in_lineage"):
            try:
                sc = variantcall.site_calls(run, r["pos"], r["alt"], max_n)
                r["callers"] = {t: f"{v['called']}/{v['judged']}" for t, v in sc["per_caller"].items()}
                r["callers_agree"] = bool(sc["patterns"]) and all(p["all"] for p in sc["patterns"])
            except Exception:
                r["callers"], r["callers_agree"] = None, None
            if tree and (r["pos"], r["alt"]) in tree:
                t = tree[(r["pos"], r["alt"])]
                r["origins"] = {"events": t["events"], "settled": t["settled"]}
    n_by = Counter(r["class"] for r in d_rows)
    return {"country": country, "private_only": private_only, "members": len(members), "others": len(O), "same_lineage_others": len(same), "lineages": dict(glin), "profile": per,
            "shared_by_all": len(shared), "classes": dict(n_by), "distinct": d_rows[:top], "n_distinct_total": len(d_rows), "lost": lost[:20], "nearest": near, "nearest_ties": ties, "countries": _countries(meta),
            "no_test": len(members) < MIN_FOR_TEST, "reference": mutprofile.summary(run).get("reference", "")}


def _countries(meta: dict) -> list[str]:
    c = Counter((geo.resolve(m["country"]) or {}).get("name") or m["country"] for m in meta.values() if m["country"])
    return [k for k, _n in c.most_common(60)]
