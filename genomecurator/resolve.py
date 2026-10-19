"""Resolving what the sources disagree about, one field at a time.

The merge picks a value by source priority and says nothing about the others. Here the person sees, for one genome, every field on which
the sources differ (what each source says, and whether it is only written differently or really different), gets a *suggestion* with
the reason and a confidence, and chooses. A choice is stored (`field_choice`), applied to the open changeset items at once, and honoured
by every later check (`curate.merge`), so a decision is made once.

The suggestion rules, in the order they apply:
  1. a value that cannot be right is left out (an impossible date: before 1900 or in the future)
  2. length / version / sequence hash: NCBI's value, because the sequence is downloaded from NCBI and its length verified against it
  3. values that are the same fact written differently form one group; the best-written one is suggested (the most specific date, the
     fullest place, the scientific host name, the standard country name)            -> confidence HIGH
  4. values that really differ: the group most sources agree on                        -> confidence MEDIUM (e.g. 2 of 3)
  5. no majority: the source priority decides, and the person is told to look          -> confidence LOW"""
import json
import re
from datetime import date, datetime
from itertools import combinations

from . import config, reconcile, seqcheck
from .db import audit

LABEL = {"ncbi": "NCBI", "bvbrc": "BV-BRC", "pathoplexus": "Pathoplexus", "nextclade": "Nextclade", "custom": "직접 입력", "manual": "직접 입력"}
CONF_RANK = {"high": 3, "medium": 2, "low": 1}
CONF_MIN = {"high": 3, "medium": 2, "low": 1}
NUMERIC = {"accession_version", "genome_length"}


def _s(v) -> str:
    return "; ".join(str(x) for x in v) if isinstance(v, (list, tuple, set)) else str(v).strip()


def _prio(sources: list[str]) -> int:
    order = config.VALUE_SOURCES
    return min((order.index(s) for s in sources if s in order), default=99)


def _date_ok(v: str) -> bool:
    p = reconcile._date_parts(v)
    if not p:
        return False
    year = int(p[0])
    if not 1900 <= year <= date.today().year + 1:
        return False
    return p[1] is None or 1 <= int(p[1]) <= 12


def _date_detail(v: str) -> tuple:
    """How informative a date is: parts present, then written in the ISO standard, then length."""
    p = reconcile._date_parts(v)
    return (sum(1 for x in p[:3] if x) if p else 0, bool(re.fullmatch(r"\d{4}(-\d{2}){0,2}", v.strip())), len(v))


def _is_binomial(v: str) -> bool:
    return bool(re.fullmatch(r"[A-Z][a-z]+ [a-z]+(?: [a-z]+)?", v.strip()))


def _canonical(kind: str, cluster: list[dict]) -> tuple[dict, str]:
    """The best-written member of a group of equivalent values, and why."""
    if len(cluster) == 1:
        return cluster[0], ""
    if kind == "date":
        return max(cluster, key=lambda g: (_date_detail(g["value"]), -_prio(g["sources"]))), "가장 구체적인 날짜"
    if kind in ("place", "contain"):
        return max(cluster, key=lambda g: (len(g["value"]), -_prio(g["sources"]))), "더 구체적인(자세한) 표기"
    if kind == "host":
        sci = [g for g in cluster if _is_binomial(g["value"].split(";")[0])]
        if sci:
            return min(sci, key=lambda g: _prio(g["sources"])), "학명(속·종) 표기"
    if kind == "country":
        std = [g for g in cluster if g["value"].casefold() not in reconcile.COUNTRY_ALIAS]
        if std:
            return min(std, key=lambda g: _prio(g["sources"])), "표준 국가명"
    return min(cluster, key=lambda g: _prio(g["sources"])), "소스 우선순위"


def _clusters(kind: str, groups: list[dict]) -> list[list[dict]]:
    """Groups of values that are the same fact (equal or only written differently)."""
    parent = list(range(len(groups)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i, j in combinations(range(len(groups)), 2):
        if reconcile.classify(kind, groups[i]["value"], groups[j]["value"]) != "conflict":
            parent[find(i)] = find(j)
    out: dict[int, list[dict]] = {}
    for i, g in enumerate(groups):
        out.setdefault(find(i), []).append(g)
    return list(out.values())


def suggest(field: str, groups: list[dict], seq: dict | None = None) -> dict:
    """-> {value, sources, confidence, reason}; `groups` = [{value, sources}] with distinct values."""
    kind = reconcile.FIELDS.get(field, "text")
    seq = seq or {}
    cands, note = groups, ""
    if kind == "date":
        ok = [g for g in groups if _date_ok(g["value"])]
        if ok and len(ok) < len(groups):
            cands, note = ok, " (불가능한 날짜 값은 제외)"
    if field in seqcheck.SEQ_FIELDS:
        g = next((g for g in cands if "ncbi" in g["sources"]), min(cands, key=lambda g: _prio(g["sources"])))
        same = [s for s, v in seq.items() if v == "same" and s != "ncbi"]
        why = "NCBI 값입니다. 서열을 NCBI에서 받아 길이를 검증하므로 NCBI 기준이 가장 안전합니다"
        if same:
            why += " (" + ", ".join(LABEL[s] for s in same) + "는 서열이 NCBI와 같아 표기·버전 지연일 뿐입니다)"
        return {**g, "confidence": "high", "reason": why}
    nc = next((g for g in cands if "nextclade" in g["sources"]), None)
    if field in ("clade", "lineage") and nc is not None and len(cands) > 1:
        others = [s for g in cands for s in g["sources"] if s != "nextclade" and g is not nc]
        agree = [s for s in nc["sources"] if s != "nextclade"]
        why = "Nextclade가 저장된 서열을 데이터셋 기준으로 분류한 값입니다 (서열에서 계산하므로 소스의 표기보다 근거가 확실합니다)"
        if agree:
            return {**nc, "confidence": "high", "reason": why + f" — {', '.join(LABEL[s] for s in agree)}도 같은 값입니다"}
        return {**nc, "confidence": "medium", "reason": why + f". {', '.join(LABEL[s] for s in others)}와 다르니 의심되면 소스 값을 고르세요"}
    clusters = _clusters(kind, cands)
    weight = lambda cl: sum(len(g["sources"]) for g in cl)
    top = max(weight(c) for c in clusters)
    leaders = [c for c in clusters if weight(c) == top]
    if len(clusters) == 1:
        g, how = _canonical(kind, clusters[0])
        norm = {re.sub(r"[\W_]+", "", x["value"].casefold()) for x in clusters[0]}
        # only punctuation / case / spacing differs -> certain; one value merely contains the other (more detail) -> worth a glance
        sure = len(norm) == 1 or kind in ("date", "country", "host")
        return {**g, "confidence": "high" if sure else "medium",
                "reason": ("같은 사실을 다르게 적은 것입니다" if sure else "한쪽이 다른 쪽을 포함합니다(더 구체적인 쪽을 제안)") + (f" — {how}을(를) 제안합니다" if how and sure else "") + note}
    if len(leaders) == 1 and top >= 2:
        g, how = _canonical(kind, leaders[0])
        n_all = sum(len(x["sources"]) for x in cands)
        return {**g, "confidence": "medium", "reason": f"{n_all}개 소스 중 {top}개가 같은 값입니다" + (f" ({how})" if how and len(leaders[0]) > 1 else "") + note}
    best = min(cands, key=lambda g: _prio(g["sources"]))
    pri = " > ".join(LABEL[s] for s in config.SOURCE_PRIORITY)
    return {**best, "confidence": "low", "reason": f"소스끼리 모두 달라 우선순위({pri})로 제안합니다. 직접 확인해서 고르세요" + note}


def analyse(per_src: dict[str, dict], current: dict, prov: dict, chosen: dict, seq: dict | None = None) -> list[dict]:
    """Every compared field on which the sources really differ for one accession.
    per_src {source: record}, current = the record that would be applied, chosen {field: {source, value}} saved choices."""
    out = []
    for field, kind in reconcile.FIELDS.items():
        vals = {s: _s(rec[field]) for s, rec in per_src.items() if rec.get(field) not in (None, "", [])}
        if len(set(vals.values())) < 2:
            continue
        verdicts = [reconcile.classify(kind, vals[a], vals[b]) for a, b in combinations(vals, 2)]
        if all(v == "match" for v in verdicts):
            continue
        worst = "conflict" if "conflict" in verdicts else "cosmetic"
        groups: dict[str, list[str]] = {}
        for s in config.VALUE_SOURCES:
            if s in vals:
                groups.setdefault(vals[s], []).append(s)
        glist = [{"value": v, "sources": ss} for v, ss in groups.items()]
        sug = suggest(field, glist, {s: v for s, v in (seq or {}).items()})
        out.append({"field": field, "kind": worst, "groups": glist, "suggestion": sug,
                    "current": {"value": current.get(field), "source": prov.get(field)}, "chosen": chosen.get(field),
                    "seq": {s: v for s, v in (seq or {}).items() if field in seqcheck.SEQ_FIELDS}})
    out.sort(key=lambda f: (f["kind"] != "conflict", list(reconcile.FIELDS).index(f["field"])))
    return out


# ------------------------------------------------------------------------------------------------ persistence
def load_choices(con, watch_id: int) -> dict[str, dict[str, dict]]:
    """{accession: {field: {source, value}}}"""
    out: dict[str, dict[str, dict]] = {}
    for acc, field, src, val in con.execute("SELECT accession, field, source, value FROM field_choice WHERE watch_id=?", (watch_id,)):
        out.setdefault(acc, {})[field] = {"source": src, "value": val}
    return out


def _coerce(field: str, value: str):
    if field in NUMERIC:
        try:
            return int(float(value))
        except ValueError:
            return value
    return value


def default_value(per_src: dict[str, dict], field: str):
    """What the plain priority merge gives for this field -> (value, source) or (None, None)."""
    for s in config.VALUE_SOURCES:
        v = per_src.get(s, {}).get(field)
        if isinstance(v, (list, tuple, set)):
            v = "; ".join(str(x) for x in v)
        if v not in (None, ""):
            return v, s
    return None, None


def _source_records(con, watch_id: int, accession: str) -> dict[str, dict]:
    return {s: json.loads(p) for s, p in con.execute("SELECT source, payload_json FROM source_record WHERE watch_id=? AND accession=?", (watch_id, accession))}


def save(con, watch_id: int, accession: str, choices: list[dict], auto: bool = False) -> dict:
    """Store choices ({field, source, value}; source None = back to the default) and apply them to the open changeset items.
    The caller commits. Returns {saved, reset, items}."""
    from . import curate
    per_src = _source_records(con, watch_id, accession)
    saved = reset = 0
    plan: dict[str, tuple] = {}   # field -> (value, source label for provenance)
    for ch in choices:
        field, src = ch.get("field"), ch.get("source")
        if field not in reconcile.FIELDS:
            raise ValueError(f"선택할 수 없는 필드입니다: {field}")
        if not src:   # back to the priority merge
            con.execute("DELETE FROM field_choice WHERE watch_id=? AND accession=? AND field=?", (watch_id, accession, field))
            v, s = default_value(per_src, field)
            plan[field] = (v, s)
            reset += 1
            continue
        if src == "custom":
            val = _s(ch.get("value", ""))
            if not val:
                raise ValueError("직접 입력한 값이 비어 있습니다.")
            plan[field] = (_coerce(field, val), "manual")
        elif src in per_src and per_src[src].get(field) not in (None, ""):
            val = _s(per_src[src][field])
            plan[field] = (per_src[src][field] if not isinstance(per_src[src][field], (list, tuple, set)) else val, src)
        else:
            raise ValueError(f"{LABEL.get(src, src)}에는 {field} 값이 없습니다.")
        con.execute("INSERT OR REPLACE INTO field_choice(watch_id, accession, field, source, value, auto, decided_at) VALUES(?,?,?,?,?,?,?)",
                    (watch_id, accession, field, src, _s(plan[field][0]), 1 if auto else 0, datetime.now().isoformat(timespec="seconds")))
        saved += 1
    items = 0
    for it in con.execute("SELECT i.* FROM changeset_item i JOIN changeset c ON c.id=i.changeset_id WHERE c.watch_id=? AND i.accession=? "
                          "AND c.status IN ('PENDING','READY','BLOCKED') AND i.decision!='applied' AND i.kind IN ('new','changed')", (watch_id, accession)).fetchall():
        rec, prov = json.loads(it["new_record_json"] or "{}"), json.loads(it["provenance_json"] or "{}")
        for field, (v, s) in plan.items():
            if v in (None, ""):
                rec.pop(field, None)
                prov.pop(field, None)
            else:
                rec[field], prov[field] = v, s
        old = con.execute("SELECT * FROM genome WHERE accession=?", (accession,)).fetchone()
        kind = it["kind"]
        diff = curate._diff(dict(old), rec) if old else {k: [None, rec.get(k)] for k in curate.COMPARE if rec.get(k) is not None}
        con.execute("UPDATE changeset_item SET new_record_json=?, provenance_json=?, diff_json=?, fingerprint=? WHERE id=?",
                    (json.dumps(rec), json.dumps(prov), json.dumps(diff), curate._fingerprint(kind, accession, rec), it["id"]))
        items += 1
    audit(con, "resolve", f"{accession}: {saved}개 선택, {reset}개 기본값 복귀" + (" (추천 일괄)" if auto else ""))
    return {"saved": saved, "reset": reset, "items": items}


def item_view(con, changeset_id: int, item_id: int) -> dict:
    """Everything the review panel needs for one changeset item."""
    it = con.execute("SELECT i.*, c.watch_id FROM changeset_item i JOIN changeset c ON c.id=i.changeset_id WHERE i.id=? AND c.id=?", (item_id, changeset_id)).fetchone()
    if not it:
        raise ValueError("항목을 찾을 수 없습니다.")
    acc, wid = it["accession"], it["watch_id"]
    per_src = _source_records(con, wid, acc)
    rec, prov = json.loads(it["new_record_json"] or "{}"), json.loads(it["provenance_json"] or "{}")
    chosen = load_choices(con, wid).get(acc, {})
    seq = seqcheck.verdicts(wid).get(acc, {}).get("sources", {})
    fields = analyse(per_src, rec, prov, chosen, seq)      # chosen fields stay in the list so a choice can be changed or undone
    return {"accession": acc, "item": item_id, "fields": fields, "sources": {s: LABEL[s] for s in per_src}, "have_sources": sorted(per_src)}


def preview_auto(con, changeset_id: int, minimum: str = "medium", item_ids: list[int] | None = None) -> dict:
    """What 'apply the suggestions' would do for the whole changeset, or only for item_ids (the rows the person selected):
    counts by confidence and the items to apply (not yet chosen only)."""
    cs = con.execute("SELECT watch_id FROM changeset WHERE id=?", (changeset_id,)).fetchone()
    if not cs:
        raise ValueError("changeset을 찾을 수 없습니다.")
    wid = cs["watch_id"]
    choices = load_choices(con, wid)
    seqv = seqcheck.verdicts(wid)
    per_acc: dict[str, list[dict]] = {}
    counts = {"high": 0, "medium": 0, "low": 0}
    items = 0
    only = {int(i) for i in item_ids} if item_ids is not None else None
    for it in con.execute("SELECT id, accession, new_record_json, provenance_json FROM changeset_item WHERE changeset_id=? AND kind IN ('new','changed') AND decision!='applied'", (changeset_id,)):
        if only is not None and it["id"] not in only:
            continue
        acc = it["accession"]
        per_src = _source_records(con, wid, acc)
        if len(per_src) < 2:
            continue
        rec, prov = json.loads(it["new_record_json"] or "{}"), json.loads(it["provenance_json"] or "{}")
        fs = [f for f in analyse(per_src, rec, prov, choices.get(acc, {}), seqv.get(acc, {}).get("sources", {})) if f["field"] not in choices.get(acc, {})]
        if not fs:
            continue
        items += 1
        for f in fs:
            counts[f["suggestion"]["confidence"]] += 1
            if CONF_RANK[f["suggestion"]["confidence"]] >= CONF_RANK[minimum]:
                same = _s(rec.get(f["field"], "")) == f["suggestion"]["value"]
                # recorded even when it equals the current value, so the item counts as decided
                per_acc.setdefault(acc, []).append({"field": f["field"], "source": f["suggestion"]["sources"][0], "value": f["suggestion"]["value"], "same": same})
    return {"watch_id": wid, "items": items, "counts": counts, "accessions": len(per_acc), "minimum": minimum, "per_acc": per_acc,
            "fields": sum(len(v) for v in per_acc.values()), "will_change": sum(1 for v in per_acc.values() for c in v if not c["same"])}


def apply_auto(con, changeset_id: int, minimum: str = "medium", item_ids: list[int] | None = None) -> dict:
    pv = preview_auto(con, changeset_id, minimum, item_ids)
    done = 0
    for acc, chs in pv["per_acc"].items():
        save(con, pv["watch_id"], acc, [{k: c[k] for k in ("field", "source", "value")} for c in chs], auto=True)
        done += len(chs)
    return {**{k: pv[k] for k in ("items", "counts", "accessions", "will_change")}, "applied": done}
