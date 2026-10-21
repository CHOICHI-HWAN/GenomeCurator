"""A Gantt chart of the mutations: for every mutation, when it was observed - from the collection dates of the genomes that carry it.

One row per mutation (position + new base), the x axis is time in ISO weeks or months, a cell is coloured by how many genomes carrying the mutation were collected in that period, and the rows are
ordered by **when each mutation was first observed** (the earliest first), so the picture reads from the mutations that appeared first at the top to the latest at the bottom. The same filters as the
genome-by-genome picture apply (lineage, clade, gene or region, APOBEC3-type only, drop what nearly every genome carries).

What "first observed" means: the earliest collection date among the sampled genomes that carry the mutation - not when the mutation arose. Sampling is uneven in time and place, so a mutation that arose
earlier but was sequenced later appears later, and a mutation present in the first genomes sampled has "first observed" at the start of the sampling period whatever its real age. A genome with only a
month for a date counts in the monthly view only (a week cannot be put on it)."""
from collections import Counter, defaultdict

import logging

from . import continents, insight, lineages, mutprofile, origins, outbreak, trend

_log = logging.getLogger("genomecurator.mutgantt")


def _key_date(kind: str, shown: str, month: str) -> tuple:
    """Sort key of a collection date: a day sorts by its date, a month-only date as the first of its month."""
    return (month, shown if kind == "day" else month + "-00")


_CONT: dict[str, str] = {}


def _cont(country: str) -> str:
    """The continent of a country name as written in the metadata ('' when unknown); cached."""
    if country not in _CONT:
        _CONT[country] = continents.group_of(country, "continent")
    return _CONT[country]


def defining_map(run: dict, meta: dict, fa) -> dict[tuple, list[str]]:
    """{(position in the run's reference, new base): [lineages that this mutation defines]} from the lineage designation (the mpxv-lineages definitions), for the lineages of the run's genomes.
    The lineage listed is the one that *introduced* the mutation (a sub-lineage inherits it and is not listed again). Empty when the definitions are not downloaded."""
    out: dict[tuple, set] = {}
    try:
        if not (lineages.status()["ok"] and fa):
            return {}
        to_run = lineages.defining_in_run_coordinates(run, fa)
        for ob, lin in {(m.get("outbreak", ""), m["lineage"]) for m in meta.values() if m["lineage"] and m.get("outbreak")}:
            for pos, nuc, by in lineages.defining_snps(ob, lin):
                p = to_run.get(pos) if to_run is not None else pos
                if p:
                    out.setdefault((p, nuc), set()).add(by)
    except Exception:
        _log.warning("lineage definitions could not be placed on the reference", exc_info=True)
        return {}
    return {k: sorted(v) for k, v in out.items()}


RECENT_MONTHS = 6
FIXED_SHARE = 0.9


def _definitional(meta: dict, s: str, introducers: list[str]) -> bool:
    """Does genome `s` carry the mutation because its lineage (or an ancestor lineage) is defined by it?"""
    m = meta[s]
    if not introducers or not m["lineage"]:
        return False
    names = {c["name"] for c in lineages.chain(m.get("outbreak", ""), m["lineage"])}
    return bool(names & set(introducers))


def adaptive_traits(k, cs, ds_months, pool_by_month, recent, meta, dmap, row, n_pool=0, orig=None, have_tree=False, tr=None) -> dict:
    """Exploratory evidence that a mutation is *kept* by the virus, from the sampled genomes: it is seen in many months and still in the most recent ones, it was rare early and common now (a spread rather than an ancestral change), it turns up in several lineages that it does not
    define (an independent origin is likely) and in several countries, it changes the protein, and its share among the genomes sampled stays up. Points are added up; this ranks candidates, it does not
    show selection - a mutation that rides on a successful branch (hitchhiking) or sits in a heavily sampled place scores high as well, and the mutation's real age and origin are not known."""
    months = Counter(ds_months)
    obs = len(months)
    last_seen = max(months) if months else ""
    persistent = obs >= 3 and last_seen in recent
    rec_n = sum(v for m, v in months.items() if m in recent)
    rec_tot = sum(pool_by_month.get(m, 0) for m in recent)
    fr = [months[m] / pool_by_month[m] for m in months if pool_by_month.get(m)]
    peak = max(fr) if fr else 0.0
    rec_freq = rec_n / rec_tot if rec_tot else 0.0
    sustained = rec_n >= 3 and peak > 0 and rec_freq >= 0.5 * peak
    allm = sorted(pool_by_month)
    early = set(allm[:max(3, len(allm) // 3)])
    e_tot = sum(pool_by_month[m] for m in early)
    early_freq = sum(v for m, v in months.items() if m in early) / e_tot if e_tot else 0.0
    rising = early_freq < 0.5 and rec_freq >= 0.5 and rec_n >= 3
    intro = dmap.get(k, [])
    nondef = Counter(meta[s]["lineage"] or "(없음)" for s in cs if not _definitional(meta, s, intro))
    n_lin = sum(1 for v in nondef.values() if v >= 2) if dmap else None
    n_ctry = len({meta[s]["country"] for s in cs if meta[s]["country"]})
    syn = (row["effect"] or "").lower().startswith("syn")
    changes = bool(row["gene"]) and not syn and row["effect"] not in ("", "intergenic")
    score = 0
    why = []
    o = (orig or {"events": 0, "settled": 0, "singles": 0, "tips": 0}) if have_tree else None
    if persistent:
        if have_tree:                                                 # independent origins on the tree (docs/ADAPTIVE.md E7); the lineage count below is then only information
            if o["settled"] >= 2:
                score += 2 + (1 if o["settled"] >= 3 else 0)
                why.append(f"트리에서 독립적으로 {o['settled']}번 생김(각각 자손 2건 이상)")
                qn = o.get("q_neutral")
                if qn is not None and qn < 0.05:
                    score += 1
                    why.append(f"같은 종류의 중립 변이보다 독립 출현이 많음(q={qn:.3f})")
                elif qn is not None:
                    why.append(f"독립 출현 횟수는 중립 변이에서도 흔한 수준(q={qn:.2f})")
            elif o["settled"] == 1:
                why.append("트리에서 한 번만 생김(한 가지가 퍼진 것일 수 있음)")
            elif o["events"] == 0:
                why.append("트리 가지에 없음(트리 뿌리 이전에 생김 — 기원은 이 트리 밖)")
        elif n_lin is not None and n_lin >= 2:
            score += 2 + (1 if n_lin >= 3 else 0)
            why.append(f"정의하지 않는 lineage {n_lin}곳에서 나타남(트리가 없어 lineage 이름으로 추정)")
        if n_ctry >= 3:
            score += 1
            why.append(f"{n_ctry}개 나라")
        if changes:
            score += 1
            why.append("단백질을 바꿈")
        if sustained:
            score += 1
            why.append("최근에도 비율 유지")
        if rising:
            score += 2
            why.append(f"초기 {early_freq:.0%} → 최근 {rec_freq:.0%}로 퍼짐")
        if tr and tr["lo"] > 0:
            score += 1
            why.append(f"나라 구성을 고정해도 증가(월 로그오즈 +{tr['slope']:.2f}, 95% {tr['lo']:+.2f}~{tr['hi']:+.2f})")
    defining = bool(intro)
    fixed = n_pool >= 20 and (len(cs) / n_pool >= FIXED_SHARE or (early_freq >= 0.8 and rec_freq >= 0.8))                       # carried by nearly every genome of the set: ancestral to it (a clade marker against the reference), not a change within it
    tier = "fixed" if fixed else "defining" if defining else ("strong" if score >= 6 and rising else "candidate" if score >= 4 else "")
    if not tier and have_tree and o["settled"] >= 2:
        tier = "recurrent"                                           # arose independently on the tree at least twice, whatever its persistence
    if fixed:
        score, why = 0, ["거의 모든 유전체가 가짐(이 데이터셋의 조상형)"]
    if tier in ("strong", "candidate", "recurrent") and row.get("profile_of_interest", "") not in ("", "other"):
        tier = "apobec"                                               # APOBEC3-type changes are edits by the host: lasting and widespread, but no evidence of selection
    evidence = {"persist": persistent, "recur": bool(have_tree and o["settled"] >= 2) if have_tree else bool(n_lin is not None and n_lin >= 2), "neutral": bool(have_tree and (o.get("q_neutral") is not None and o["q_neutral"] < 0.05)),
                "spread": n_ctry >= 3, "protein": changes, "rise": rising, "adjusted": bool(tr and tr["lo"] > 0)}
    return {"score": score, "tier": tier, "evidence": evidence, "months_seen": obs, "last_seen": last_seen, "persistent": persistent, "lineages": n_lin, "origins": ({**o, "source": "tree"} if have_tree else {"source": "lineage"}), "countries": n_ctry, "recent_n": rec_n, "peak_freq": round(peak, 4), "early_freq": round(early_freq, 4), "trend": ({k_: round(v_, 4) for k_, v_ in tr.items()} if tr else None), "rising": rising,
            "recent_freq": round(rec_freq, 4), "changes_protein": changes, "why": why}



def co_travel(rows: list[dict], min_j: float = 0.6) -> list[dict]:
    """Pairs of the given mutations that sit on (nearly) the same genomes: Jaccard similarity of their carrier sets, `a` the one first observed earlier. Mutations that always travel together were most
    likely gained together (one event or one branch), so they are one piece of evidence, not several."""
    out = []
    for i, x in enumerate(rows):
        for y in rows[i + 1:]:
            inter = len(x["_cs"] & y["_cs"])
            if not inter:
                continue
            j = inter / len(x["_cs"] | y["_cs"])
            if j >= min_j:
                a, b = sorted((x, y), key=lambda r: (r["first"], r["pos"]))
                out.append({"a": f"{a['pos']}{a['alt']}", "b": f"{b['pos']}{b['alt']}", "j": round(j, 3), "shared": inter})
    return out

def build(run: dict, bin_: str = "month", lineage_set=(), clade_set=(), gene: str = "", start: int = 0, end: int = 0, apobec_only: bool = False, hide_common: float = 0.0,
          min_carriers: int = 3, max_rows: int = 150, sort: str = "first", max_n: float | None = 0.05, defining: str = "all") -> dict:
    sm = mutprofile.summary(run)
    ref = sm.get("reference") or ""
    fa, gff = mutprofile.reference_by_name(ref) if ref else (None, None)
    refname, refseq = mutprofile.read_fasta(fa) if fa else ("", "")
    if gene:
        genes = {c["gene"]: c for c in mutprofile.load_cds(gff)}
        if gene not in genes:
            raise ValueError(f"유전자 {gene}을(를) 참조 주석에서 찾을 수 없습니다.")
        ex = genes[gene]["exons"]
        start, end = max(1, min(a for a, _ in ex) - 200), max(b for _, b in ex) + 200
    start, end = start or 1, end or len(refseq) or 10 ** 9
    bin_ = bin_ if bin_ in ("week", "month") else "month"
    meta = mutprofile._run_meta(run)
    dmap = defining_map(run, meta, fa)
    lset, cset = set(lineage_set), set(clade_set)
    rows = insight._mutations(run, max_n)                                     # ITR left out, genomes with too many N's left out
    carriers: dict[tuple, set] = defaultdict(set)
    info: dict[tuple, dict] = {}
    pool: set[str] = set()
    for r in rows:
        s = r["sample"]
        m = meta.get(s)
        if m is None or (lset and m["lineage"] not in lset) or (cset and m["clade"] not in cset):
            continue
        pool.add(s)
        pos = int(r["pos"])
        if not (start <= pos <= end) or (apobec_only and r["profile_of_interest"] in ("", "other")):
            continue
        k = (pos, r["alt"])
        carriers[k].add(s)
        info.setdefault(k, r)
    dated: dict[str, tuple] = {}
    for s in pool:
        kind, week, month, shown = lineages.parse_date(meta[s]["date"])
        if kind:
            dated[s] = (kind, week, month, shown)
    n_pool = max(1, len(pool))
    pool_by_month = Counter(v[2] for v in dated.values())
    recent = set(sorted(pool_by_month)[-RECENT_MONTHS:]) if pool_by_month else set()
    cats = outbreak.categories()
    origs = origins.independent_origins(run)
    month_of = {s_: {m: i for i, m in enumerate(sorted(pool_by_month))}[dated[s_][2]] for s_ in dated}
    tctx = trend.context(sorted(dated), month_of, {s_: meta[s_]["country"] for s_ in dated})
    out = []
    for k, cs in carriers.items():
        if len(cs) < min_carriers or (hide_common and len(cs) / n_pool >= hide_common):
            continue
        if (defining == "only" and k not in dmap) or (defining == "none" and k in dmap):
            continue
        ds = [(dated[s], s) for s in cs if s in dated]
        if not ds:
            continue
        ds.sort(key=lambda x: _key_date(x[0][0], x[0][3], x[0][2]))
        cells: Counter = Counter()
        for (kind, week, month, _shown), _s in ds:
            if bin_ == "month":
                cells[month] += 1
            elif kind == "day":
                cells[week] += 1
        r = info[k]
        first, last = ds[0][0], ds[-1][0]
        out.append({"key": f"{k[0]}{k[1]}", "label": f"{r['ref']}{k[0]}{k[1]}", "pos": k[0], "ref": r["ref"], "alt": k[1], "gene": r["gene"], "aa": r["aa_change"], "effect": r["effect"],
                    "apobec": int(r["profile_of_interest"] not in ("", "other")), "defines": dmap.get(k, []), "n": len(cs), "n_dated": len(ds), "first": first[3], "last": last[3], "first_month": first[2], "last_month": last[2],
                    "first_week": first[1], "last_week": last[1], "first_lineage": meta[ds[0][1]]["lineage"], "first_country": meta[ds[0][1]]["country"], "cells": dict(cells),
                    "countries": dict(Counter(meta[s]["country"] for _d, s in ds if meta[s]["country"]).most_common(3)),
                    "by_lineage": dict(Counter(meta[s]["lineage"] or "(없음)" for s in cs)), "by_country": dict(Counter(meta[s]["country"] for s in cs if meta[s]["country"])),
                    "by_continent": dict(Counter(_cont(meta[s]["country"]) for s in cs if _cont(meta[s]["country"]))),
                    "_cs": cs,
                    "category": (cats.get(r["gene"], {}).get("category") or ("(유전자 밖)" if not r["gene"] else "(분류 없음)")), "function": cats.get(r["gene"], {}).get("function", ""),
                    "adaptive": adaptive_traits(k, cs, [d[2] for d, _s in ds], pool_by_month, recent, meta, dmap, r, len(pool), origs.get((k[0], k[1])) if origs else None, origs is not None, trend.carrier_trend(cs, tctx)),
                    "_key": _key_date(first[0], first[3], first[2])})
    total_rows = len(out)
    if sort == "pos":
        out.sort(key=lambda x: (x["pos"], x["alt"]))
    elif sort == "n":
        out.sort(key=lambda x: (-x["n"], x["_key"], x["pos"]))
    else:
        out.sort(key=lambda x: (x["_key"], x["pos"], x["alt"]))
    out = out[:max_rows]
    links = co_travel([x for x in out if x["adaptive"]["tier"] in ("strong", "candidate", "recurrent", "apobec")])
    for x in out:
        x.pop("_key")
        x.pop("_cs")
    all_dates = [dated[s] for s in dated]
    if bin_ == "week":
        wk = sorted({d[1] for d in all_dates if d[0] == "day"})
        periods = lineages._weeks(wk[0], wk[-1]) if wk else []
        starts = [lineages.week_monday(k).isoformat() for k in periods]
    else:
        mo = sorted({d[2] for d in all_dates})
        periods = lineages._months(mo[0], mo[-1]) if mo else []
        starts = [k + "-01" for k in periods]
    return {"reference": refname or ref, "bin": bin_, "periods": periods, "starts": starts, "rows": out, "total_rows": total_rows, "shown_rows": len(out), "genomes": len(pool), "dated": len(dated),
            "month_only": sum(1 for d in all_dates if d[0] == "month"), "no_date": len(pool) - len(dated), "start": start, "end": end, "sort": sort, "min_carriers": min_carriers, "max_n": max_n, "definitions": bool(dmap), "origins_from_tree": origs is not None, "defining": defining,
            "n_defining_shown": sum(1 for x in out if x["defines"]), "adaptive_links": links}


def series(run: dict, pos: int, alt: str, bin_: str = "month", max_n: float | None = 0.05) -> dict:
    """Per month or ISO week: how many of the sampled genomes carry one mutation (stacked by lineage), and how many genomes were sampled in that period - so a rise can be told from heavier sampling.
    Genomes with only a month for a date appear in the monthly view only; undated genomes are counted separately."""
    bin_ = bin_ if bin_ in ("week", "month") else "month"
    meta = mutprofile._run_meta(run)
    carriers = {r["sample"] for r in insight._mutations(run, max_n) if int(r["pos"]) == pos and r["alt"] == alt}
    sampled: Counter = Counter()
    by_lin: dict[str, Counter] = defaultdict(Counter)
    undated = 0
    for s, m in meta.items():
        kind, week, month, _shown = lineages.parse_date(m["date"])
        key = month if bin_ == "month" else (week if kind == "day" else None)
        if not kind:
            undated += 1
            continue
        if key is None:
            continue
        sampled[key] += 1
        if s in carriers:
            by_lin[m["lineage"] or "(none)"][key] += 1
    ckeys = sorted({k for c in by_lin.values() for k in c})                                                       # the span of the carriers: the same span the spread timeline draws, so the two line up
    keys = [k for k in sorted(sampled) if ckeys and ckeys[0] <= k <= ckeys[-1]] or ckeys
    if bin_ == "week":
        periods = lineages._weeks(keys[0], keys[-1]) if keys else []
        starts = [lineages.week_monday(k).isoformat() for k in periods]
    else:
        periods = lineages._months(keys[0], keys[-1]) if keys else []
        starts = [k + "-01" for k in periods]
    tot = {lin: sum(c.values()) for lin, c in by_lin.items()}
    order = sorted(by_lin, key=lambda l: -tot[l])
    return {"bin": bin_, "periods": periods, "starts": starts, "sampled": [sampled.get(p, 0) for p in periods],
            "lineages": [{"lineage": l, "counts": [by_lin[l].get(p, 0) for p in periods], "n": tot[l]} for l in order],
            "carriers": [sum(by_lin[l].get(p, 0) for l in by_lin) for p in periods], "n_carriers": len(carriers), "undated": undated}
