"""When did the genomes of a virus come into the curated table, and when were they collected?

*Intake*: the day (or ISO week, or month) a genome entered the table (`genome.date_inserted`: the moment its change set was applied; `private_genome.imported` for the private ones), plus the
changes to genomes that were already there (`genome_version.changed_at`), with the change sets that were applied marked. *Collection*: the collection dates of the same genomes, optionally only of
the genomes that came in during one intake period. The two time axes are different (the first is the history of this table, the second the history of the outbreak)."""
from collections import Counter, defaultdict

from . import lineages
from .db import connect


def _key(ts: str, bin_: str):
    """Period key of a timestamp / date text for the intake chart ('' when unusable): day -> YYYY-MM-DD, week -> ISO week, month -> YYYY-MM."""
    kind, wk, mo, shown = lineages.parse_date((ts or "")[:10])
    if not kind:
        return ""
    if bin_ == "month":
        return mo
    if bin_ == "week":
        return wk if kind == "day" else ""
    return shown if kind == "day" else ""


def _span(keys: list[str], bin_: str) -> list[str]:
    if not keys:
        return []
    if bin_ == "month":
        return lineages._months(keys[0], keys[-1])
    if bin_ == "week":
        return lineages._weeks(keys[0], keys[-1])
    from datetime import date, timedelta
    a, b = date.fromisoformat(keys[0]), date.fromisoformat(keys[-1])
    return [(a + timedelta(days=i)).isoformat() for i in range((b - a).days + 1)]


def _start(key: str, bin_: str) -> str:
    return lineages.week_monday(key).isoformat() if bin_ == "week" else (key + "-01" if bin_ == "month" else key)


def history(watch_id: int, bin_: str = "day", cbin: str = "week", sel: str = "") -> dict:
    """Intake per `bin_` (day / week / month) and the collection dates per `cbin`; `sel` = an intake period key (as in the result's `periods`) restricts the collection chart to the genomes that came in then."""
    bin_ = bin_ if bin_ in ("day", "week", "month") else "day"
    cbin = cbin if cbin in ("day", "week", "month") else "week"
    with connect() as con:
        gen = [dict(r) for r in con.execute("SELECT accession, date_inserted, collection_date FROM genome WHERE watch_id=? AND IFNULL(genome_status,'')!='Removed'", (watch_id,))]
        priv = [dict(r) for r in con.execute("SELECT id, imported, metadata_json FROM private_genome WHERE watch_id=?", (watch_id,))]
        ver = [dict(r) for r in con.execute("SELECT v.changed_at, v.changeset_id FROM genome_version v JOIN genome g ON g.genome_id=v.genome_id WHERE g.watch_id=?", (watch_id,))]
        cs = [dict(r) for r in con.execute("SELECT id, applied_at, summary_json FROM changeset WHERE watch_id=? AND status='APPLIED' AND applied_at IS NOT NULL ORDER BY applied_at", (watch_id,))]
    import json
    for p in priv:
        try:
            p["collection_date"] = (json.loads(p["metadata_json"] or "{}").get("date") or json.loads(p["metadata_json"] or "{}").get("collection_date") or "")
        except ValueError:
            p["collection_date"] = ""
    added: Counter = Counter()
    added_priv: Counter = Counter()
    changed: Counter = Counter()
    unusable = 0
    for g in gen:
        k = _key(g["date_inserted"], bin_)
        if k:
            added[k] += 1
        else:
            unusable += 1
    for p in priv:
        k = _key(p["imported"], bin_)
        if k:
            added_priv[k] += 1
    for v in ver:
        k = _key(v["changed_at"], bin_)
        if k:
            changed[k] += 1
    keys = sorted(set(added) | set(added_priv) | set(changed))
    periods = _span(keys, bin_)
    applied = []
    for c in cs:
        k = _key(c["applied_at"], bin_)
        if k:
            applied.append({"id": c["id"], "period": k, "at": c["applied_at"][:16].replace("T", " ")})
    pool = [(g["collection_date"], _key(g["date_inserted"], bin_)) for g in gen] + [(p["collection_date"], _key(p["imported"], bin_)) for p in priv]
    if sel:
        pool = [x for x in pool if x[1] == sel]
    ckeys, ccount = [], Counter()
    no_date = year_only = month_only_skipped = 0
    for d, _k in pool:
        d = (d or "").strip()
        kind, wk, mo, shown = lineages.parse_date(d)
        if not kind:
            if len(d) == 4 and d.isdigit():
                year_only += 1
            else:
                no_date += 1
            continue
        k = _key(d, cbin) if cbin != "day" else (shown if kind == "day" else "")
        if cbin == "week" and kind != "day":
            month_only_skipped += 1
            continue
        if cbin == "day" and kind != "day":
            month_only_skipped += 1
            continue
        if k:
            ccount[k] += 1
    ckeys = sorted(ccount)
    cperiods = _span(ckeys, cbin)
    return {"bin": bin_, "cbin": cbin, "periods": periods, "starts": [_start(p, bin_) for p in periods],
            "added": [added.get(p, 0) for p in periods], "added_private": [added_priv.get(p, 0) for p in periods], "changed": [changed.get(p, 0) for p in periods],
            "applied": applied, "total": len(gen) + len(priv), "private": len(priv), "unusable_intake": unusable,
            "sel": sel, "sel_n": len(pool),
            "cperiods": cperiods, "cstarts": [_start(p, cbin) for p in cperiods], "collected": [ccount.get(p, 0) for p in cperiods],
            "no_date": no_date, "year_only": year_only, "coarse_skipped": month_only_skipped}
