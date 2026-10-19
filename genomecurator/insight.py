"""Questions asked of a mutation-profile run's mutations: are APOBEC3-type mutations transmitted, which non-defining mutations are new, which mutations mark a lineage.

All three work on the consensus sequences and the SNPs of `mutations.tsv` (positions in the run's reference), and on the lineage definitions when they have been fetched.

**Are APOBEC3-type mutations transmitted?** The consensus of a genome shows the dominant virus of a sample, not the mix inside the patient. What can be told from a set of
consensus genomes is whether a mutation is *shared* - carried by more than one genome, i.e. it was passed on (at the level of the dominant virus) - or *private to one
genome* (arose in that patient or sample, or is an artefact). If APOBEC3 editing were only a within-host / sample-level event that does not get passed on, APOBEC3-type
changes would be enriched among private mutations and rare among shared ones. The share of APOBEC3-type changes is therefore compared between private and shared
mutations (by site, and - when a Nextstrain tree of the dataset exists - by branch: terminal branches carry what is private to a tip, internal branches what was inherited
by two or more). A difference is tested with a chi-square test. Limits: minor variants within a patient are invisible in a consensus (that needs reads, or several samples
from one patient); a shared mutation could also be a recurrent one at a hot spot; context is taken from the reference (sites) or the tree's root sequence (branches)."""
import json
import math
from collections import Counter, defaultdict
from datetime import date, timedelta

from . import lineages, log, mutprofile, nextstrain

_log = log.get("insight")
APO = ("GA>AA", "TC>TT")


def _mutations(run: dict, max_n: float | None = None) -> list[dict]:
    """The run's mutations without the ITR and, with `max_n`, only the genomes whose share of N / ambiguity characters is at most that."""
    rows = mutprofile._load_mutations(run)
    ok = mutprofile.good_samples(run, max_n)
    return rows if ok is None else [r for r in rows if r["sample"] in ok]


def _is_apo(row: dict) -> bool:
    return row.get("profile_of_interest", "") not in ("", "other")


def chi2_2x2(a: int, b: int, c: int, d: int) -> dict:
    """[[a, b], [c, d]] (rows: private / shared; columns: APOBEC3-type / other). chi-square (df 1, no correction), p and the odds ratio."""
    n = a + b + c + d
    if min(a + b, c + d, a + c, b + d) == 0 or n == 0:
        return {"chi2": None, "p": None, "odds_ratio": None}
    chi = n * (a * d - b * c) ** 2 / ((a + b) * (c + d) * (a + c) * (b + d))
    return {"chi2": round(chi, 2), "p": math.erfc(math.sqrt(chi / 2)), "odds_ratio": round((a * d) / (b * c), 2) if b * c else None}


def _carrier_class(k: int) -> str:
    return "1" if k == 1 else "2-3" if k <= 3 else "4-10" if k <= 10 else "11-50" if k <= 50 else ">50"


CLASSES = ["1", "2-3", "4-10", "11-50", ">50"]


def apobec_by_sites(run: dict, max_n: float | None = None) -> dict:
    """Per site (position + new base): how many genomes carry it and whether it is of the profile of interest; the share of APOBEC3-type sites / events by carrier class."""
    carriers: dict[tuple, set] = defaultdict(set)
    apo: dict[tuple, bool] = {}
    for r in _mutations(run, max_n):
        k = (int(r["pos"]), r["alt"])
        carriers[k].add(r["sample"])
        apo[k] = _is_apo(r)
    by = {c: {"sites": 0, "apobec_sites": 0, "events": 0, "apobec_events": 0} for c in CLASSES}
    for k, s in carriers.items():
        c = by[_carrier_class(len(s))]
        c["sites"] += 1
        c["events"] += len(s)
        if apo[k]:
            c["apobec_sites"] += 1
            c["apobec_events"] += len(s)
    for c in by.values():
        c["apobec_share_sites"] = round(c["apobec_sites"] / c["sites"], 4) if c["sites"] else None
    priv = by["1"]
    shared = {x: sum(by[c][x] for c in CLASSES[1:-1]) for x in ("sites", "apobec_sites")}      # carried by 2..50 genomes: passed on within the dataset
    fixed = by[CLASSES[-1]]                                                                  # >50: differences from the reference that (nearly) every genome has - not events of this dataset
    test = chi2_2x2(priv["apobec_sites"], priv["sites"] - priv["apobec_sites"], shared["apobec_sites"], shared["sites"] - shared["apobec_sites"])
    return {"classes": [{"carriers": c, **by[c]} for c in CLASSES], "private": priv, "shared": shared, "fixed": fixed, "test": test,
            "share_private": round(priv["apobec_sites"] / priv["sites"], 4) if priv["sites"] else None,
            "share_shared": round(shared["apobec_sites"] / shared["sites"], 4) if shared["sites"] else None,
            "share_fixed": fixed["apobec_share_sites"]}


def _apo_branch(ref: str, alt: str, pos: int, seq: str) -> bool:
    i = pos - 1
    if ref == "G" and alt == "A" and i + 1 < len(seq) and seq[i + 1] == "A":
        return True                                           # GA>AA
    if ref == "C" and alt == "T" and i - 1 >= 0 and seq[i - 1] == "T":
        return True                                           # TC>TT
    return False


def apobec_by_branches(run: dict) -> dict | None:
    """The same question on the Nextstrain tree of the dataset (the newest finished run): terminal branches (private to a tip) against internal branches (inherited by two or
    more genomes). Context from the tree's root sequence; the ITR ends are left out. Besides the plain comparison, three stricter variants are tested, because the plain one counts
    the same recurrent mutation many times: only distinct mutations, without mutations seen on more than one branch (homoplasy), and both. None when there is no tree."""
    import re
    runs = [r for r in nextstrain.list_runs(run["dataset_id"]) if r["tool"] == "nextstrain" and r["status"] == "done" and r["auspice_file"]]
    if not runs:
        return None
    f = nextstrain.auspice_dir() / runs[0]["auspice_file"]
    rootf = nextstrain.auspice_dir() / runs[0]["auspice_file"].replace(".json", "_root-sequence.json")
    if not f.exists() or not rootf.exists():
        return None
    tree = json.loads(f.read_text(encoding="utf-8"))["tree"]
    root = json.loads(rootf.read_text(encoding="utf-8")).get("nuc", "")
    L = len(root)
    per_branch = Counter()                                   # how many branches carry each mutation (homoplasy)
    muts: list[tuple[bool, tuple]] = []
    stack = [tree]
    while stack:
        n = stack.pop()
        kids = n.get("children") or []
        stack.extend(kids)
        for m in (n.get("branch_attrs", {}).get("mutations", {}) or {}).get("nuc", []):
            g = re.fullmatch(r"([ACGT])(\d+)([ACGT])", m)
            if not g or mutprofile.in_itr(int(g[2]), L):
                continue
            key = (g[1], int(g[2]), g[3])
            per_branch[key] += 1
            muts.append((bool(kids), key))

    def tally(distinct: bool, no_homoplasy: bool):
        seen = {True: set(), False: set()}
        c = {"internal": Counter(), "terminal": Counter()}
        for internal, key in muts:
            if no_homoplasy and per_branch[key] > 1:
                continue
            if distinct:
                if key in seen[internal]:
                    continue
                seen[internal].add(key)
            cc = c["internal" if internal else "terminal"]
            cc["n"] += 1
            cc["apobec"] += _apo_branch(key[0], key[2], key[1], root)
        t, i = c["terminal"], c["internal"]
        return {"terminal": dict(t), "internal": dict(i), "share_terminal": round(t["apobec"] / t["n"], 4) if t["n"] else None, "share_internal": round(i["apobec"] / i["n"], 4) if i["n"] else None,
                "test": chi2_2x2(t["apobec"], t["n"] - t["apobec"], i["apobec"], i["n"] - i["apobec"])}
    plain = tally(False, False)
    variants = [{"name": "전체 (같은 변이도 가지마다 셈)", **plain}, {"name": "서로 다른 변이만", **tally(True, False)}, {"name": "여러 가지에 나타난 변이(homoplasy) 제외", **tally(False, True)},
                {"name": "서로 다른 변이만 + homoplasy 제외", **tally(True, True)}]
    return {"tree": runs[0]["name"], **plain, "variants": variants, "homoplastic_mutations": sum(1 for v in per_branch.values() if v > 1), "mutations": len(per_branch)}


# ------------------------------------------------------------------------------------------------ new mutations
def _to_date(s: str):
    try:
        y, m, *rest = s.split("-")
        return date(int(y), int(m), int(rest[0]) if rest else 15)
    except (ValueError, IndexError):
        return None


def _defining_run_set(run: dict, ref_fa) -> set:
    """(position in the run's reference, new base) of every defining SNP of every lineage in the definitions (empty without definitions)."""
    if not lineages.status()["ok"]:
        return set()
    to_run = lineages.defining_in_run_coordinates(run, ref_fa)
    out = set()
    for ob, d in lineages.definitions().items():
        for lin in d:
            for pos, nuc, _by in lineages.defining_snps(ob, lin):
                p = pos if to_run is None else to_run.get(pos)
                if p:
                    out.add((p, nuc))
    return out


def emerging(run: dict, window_days: int = 120, min_carriers: int = 3, skip_apobec: bool = False, max_share: float = 0.5, max_n: float | None = None) -> dict:
    """Mutations that are not lineage-defining, carried by at least `min_carriers` genomes, not by more than `max_share` of all, ordered by first appearance. `recent` means the first
    collection date is within `window_days` of the newest date in the dataset. Mutations with the same set of carriers form a **cluster** (they travel together: a candidate
    new lineage when it spans several genomes / countries). The reference positions are the run's."""
    ref = mutprofile.summary(run).get("reference", "")
    fa, _gff = mutprofile.reference_by_name(ref) if ref else (None, None)
    gs = {g["id"]: g for g in lineages._genomes(run)}
    dates = {i: lineages.parse_date(g["date"])[3] for i, g in gs.items()}
    known = {i: _to_date(s) for i, s in dates.items() if s}
    newest = max(known.values()) if known else None
    defset = _defining_run_set(run, fa) if fa else set()
    carriers: dict[tuple, set] = defaultdict(set)
    info: dict[tuple, dict] = {}
    rows_ = _mutations(run, max_n)
    for r in rows_:
        k = (int(r["pos"]), r["alt"])
        carriers[k].add(r["sample"])
        info.setdefault(k, r)
    nall = len({r["sample"] for r in rows_}) or 1
    rows = []
    for k, cs in carriers.items():
        if len(cs) < min_carriers or len(cs) / nall > max_share or k in defset:
            continue
        r = info[k]
        if skip_apobec and _is_apo(r):
            continue
        ds_ = sorted((known[i], i) for i in cs if i in known)
        first = ds_[0][0] if ds_ else None
        rows.append({"pos": k[0], "ref": r["ref"], "alt": k[1], "gene": r["gene"], "aa": r["aa_change"], "effect": r["effect"], "apobec": _is_apo(r), "n": len(cs),
                     "countries": sorted({gs[i]["country"] for i in cs if i in gs and gs[i]["country"]}), "lineages": dict(Counter(gs[i]["lineage"] for i in cs if i in gs).most_common(4)),
                     "first": dates.get(ds_[0][1], "") if ds_ else "", "last": dates.get(ds_[-1][1], "") if ds_ else "",
                     "recent": bool(first and newest and (newest - first).days <= window_days)})
    rows.sort(key=lambda x: (not x["recent"], x["first"] or "9", -x["n"]))
    groups: dict[frozenset, list] = defaultdict(list)
    for k, cs in carriers.items():
        if min_carriers <= len(cs) and len(cs) / nall <= max_share and k not in defset:
            groups[frozenset(cs)].append(k)
    clusters = []
    for cs, ks in groups.items():
        if len(ks) < 2:
            continue
        ds_ = sorted(dates[i] for i in cs if dates.get(i))
        clusters.append({"mutations": [f"{info[k]['ref']}{k[0]}{k[1]}" for k in sorted(ks)], "genes": sorted({info[k]["gene"] for k in ks if info[k]["gene"]}), "n": len(cs),
                         "countries": sorted({gs[i]["country"] for i in cs if i in gs and gs[i]["country"]}), "lineages": dict(Counter(gs[i]["lineage"] for i in cs if i in gs).most_common(3)),
                         "first": ds_[0] if ds_ else "", "last": ds_[-1] if ds_ else "", "apobec": sum(1 for k in ks if _is_apo(info[k]))})
    clusters.sort(key=lambda c: (c["first"] or "0", c["n"]), reverse=True)
    return {"newest": newest.isoformat() if newest else "", "window_days": window_days, "min_carriers": min_carriers, "definitions": bool(defset), "reference": ref,
            "mutations": rows[:300], "n_mutations": len(rows), "clusters": clusters[:100], "n_clusters": len(clusters)}


# ------------------------------------------------------------------------------------------------ lineage signatures
def signatures(run: dict, min_n: int = 5, inside: float = 0.9, outside: float = 0.05, max_n: float | None = None) -> dict:
    """Per lineage (with its descendants counted as inside): the mutations carried by at least `inside` of its genomes and at most `outside` of all other genomes, compared with
    the official defining SNPs (how many of those are carried by the lineage's genomes), the APOBEC3-type share and the substitution types of its signature."""
    ref = mutprofile.summary(run).get("reference", "")
    fa, _g = mutprofile.reference_by_name(ref) if ref else (None, None)
    gs = {g["id"]: g for g in lineages._genomes(run)}
    by_sample: dict[str, set] = defaultdict(set)
    info: dict[tuple, dict] = {}
    for r in _mutations(run, max_n):
        k = (int(r["pos"]), r["alt"])
        by_sample[r["sample"]].add(k)
        info.setdefault(k, r)
    defs = lineages.definitions() if lineages.status()["ok"] else {}
    to_run = lineages.defining_in_run_coordinates(run, fa) if (defs and fa) else None
    try:
        nc = lineages._nc_snps(run) if defs else {}                                  # SNPs against NC_063383.1: whether a genome carries a defining SNP (the run's own SNPs miss it when the reference has the base)
        gaps = lineages.nc_unknown(run) if defs else {}                              # where a genome has no usable base (N, deletion, not covered): unknown, not absent
    except Exception:
        _log.warning("NC_063383.1 SNPs unavailable", exc_info=True)
        nc, gaps = {}, {}
    out = []
    ob_of = Counter(g["outbreak"] for g in gs.values() if g["outbreak"])
    for outbreak in sorted(ob_of):
        d = defs.get(outbreak, {})
        members: dict[str, set] = defaultdict(set)
        for i, g in gs.items():
            if g["outbreak"] == outbreak and g["lineage"] and i in by_sample:
                lin = g["lineage"]
                members[lin].add(i)
        # descendants: a lineage is inside its ancestors
        inside_of: dict[str, set] = defaultdict(set)
        for lin, ids in members.items():
            cur, seen = lin, set()
            while cur and cur not in seen:
                seen.add(cur)
                inside_of[cur] |= ids
                cur = d.get(cur, {}).get("parent", "")
        everyone = {i for ids in members.values() for i in ids}
        for lin, ids in sorted(inside_of.items()):
            if len(members.get(lin, ())) < min_n and len(ids) < min_n:
                continue
            others = everyone - ids
            cnt_in, cnt_out = Counter(), Counter()
            for i in ids:
                cnt_in.update(by_sample[i])
            for i in others:
                cnt_out.update(by_sample[i])
            sig = sorted(k for k, v in cnt_in.items() if v / len(ids) >= inside and (cnt_out[k] / len(others) if others else 0) <= outside)
            own = []
            for pos, nuc in d.get(lin, {}).get("snps", []):
                p = pos if to_run is None else (to_run.get(pos) if to_run else None)
                known = [i for i in ids if not lineages.is_unknown(gaps.get(i), pos)] if nc else list(ids)
                have = sum(1 for i in known if nc.get(i, {}).get(pos) == nuc) if nc else (cnt_in.get((p, nuc), 0) if p else None)
                own.append({"snp": f"{pos}{nuc}", "carried": (round(have / len(known), 3) if have is not None and known else None), "unknown": round(1 - len(known) / len(ids), 3) if ids else 0, "in_signature": bool(p and (p, nuc) in set(sig)), "liftable": bool(p)})
            conc = None
            if nc and d.get(lin, {}).get("snps") is not None:
                need = [(pos, nuc) for pos, nuc, _b in lineages.defining_snps(outbreak, lin)]
                mine = members.get(lin, set())
                full = [i for i in mine if need and not any(lineages.is_unknown(gaps.get(i), pos) for pos, _n in need)]
                conc = {"genomes": len(mine), "all_sites_known": len(full), "carry_all": sum(1 for i in full if all(nc.get(i, {}).get(pos) == nuc for pos, nuc in need)),
                        "with_unknown": len(mine) - len(full)} if need else None
            types = Counter(f"{info[k]['ref']}>{info[k]['alt']}" for k in sig)
            out.append({"outbreak": outbreak, "lineage": lin, "n": len(ids), "own_n": len(members.get(lin, ())), "signature_sites": len(sig),
                        "apobec_share": round(sum(1 for k in sig if _is_apo(info[k])) / len(sig), 3) if sig else None,
                        "official": own, "concordance": conc, "types": dict(types.most_common(6)),
                        "genes": dict(Counter(info[k]["gene"] or "(intergenic)" for k in sig).most_common(4)),
                        "signature": [f"{info[k]['ref']}{k[0]}{k[1]}" for k in sig[:12]]})
    return {"reference": ref, "definitions": bool(defs), "min_n": min_n, "inside": inside, "outside": outside, "lineages": out}


def insight(run: dict, **kw) -> dict:
    mn = kw.get("max_n")
    return {"apobec_sites": apobec_by_sites(run, mn), "apobec_branches": apobec_by_branches(run), "emerging": emerging(run, **{k: v for k, v in kw.items() if k in ("window_days", "min_carriers", "skip_apobec", "max_n")}),
            "signatures": signatures(run, max_n=mn), "max_n": mn, "n_filtered": _n_filtered(run, mn)}


def _n_filtered(run: dict, max_n: float | None) -> dict:
    ok = mutprofile.good_samples(run, max_n)
    allr = {r["sample"] for r in mutprofile._load_mutations(run)}
    return {"kept": len(allr & ok) if ok is not None else len(allr), "of": len(allr), "max_n": max_n}
