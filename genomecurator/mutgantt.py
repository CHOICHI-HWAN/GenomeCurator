"""A Gantt chart of the mutations: for every mutation, when it was observed - from the collection dates of the genomes that carry it.

One row per mutation (position + new base), the x axis is time in ISO weeks or months, a cell is coloured by how many genomes carrying the mutation were collected in that period, and the rows are
ordered by **when each mutation was first observed** (the earliest first), so the picture reads from the mutations that appeared first at the top to the latest at the bottom. The same filters as the
genome-by-genome picture apply (lineage, clade, gene or region, APOBEC3-type only, drop what nearly every genome carries).

What "first observed" means: the earliest collection date among the sampled genomes that carry the mutation - not when the mutation arose. Sampling is uneven in time and place, so a mutation that arose
earlier but was sequenced later appears later, and a mutation present in the first genomes sampled has "first observed" at the start of the sampling period whatever its real age. A genome with only a
month for a date counts in the monthly view only (a week cannot be put on it)."""
from collections import Counter, defaultdict

import logging

from . import insight, lineages, mutprofile

_log = logging.getLogger("genomecurator.mutgantt")


def _key_date(kind: str, shown: str, month: str) -> tuple:
    """Sort key of a collection date: a day sorts by its date, a month-only date as the first of its month."""
    return (month, shown if kind == "day" else month + "-00")


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
                    "_key": _key_date(first[0], first[3], first[2])})
    total_rows = len(out)
    if sort == "pos":
        out.sort(key=lambda x: (x["pos"], x["alt"]))
    elif sort == "n":
        out.sort(key=lambda x: (-x["n"], x["_key"], x["pos"]))
    else:
        out.sort(key=lambda x: (x["_key"], x["pos"], x["alt"]))
    out = out[:max_rows]
    for x in out:
        x.pop("_key")
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
            "month_only": sum(1 for d in all_dates if d[0] == "month"), "no_date": len(pool) - len(dated), "start": start, "end": end, "sort": sort, "min_carriers": min_carriers, "max_n": max_n, "definitions": bool(dmap), "defining": defining,
            "n_defining_shown": sum(1 for x in out if x["defines"])}


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
