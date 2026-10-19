import json
import time
from pathlib import Path

from fastapi import FastAPI, Form, Request
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse, Response, StreamingResponse
from fastapi.templating import Jinja2Templates

from . import __version__, analysis, completeness, export, qc, config, curate, jobs, log
from .db import audit_now, connect

_log = log.get("web")
_alert_cache = {"t": 0.0, "v": None}

T = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))


def _rows(sql, *a):
    con = connect()
    try:
        return [dict(r) for r in con.execute(sql, a)]
    finally:
        con.close()


def _tsv(rows: list[dict], filename: str) -> Response:
    cols = list(dict.fromkeys(k for r in rows for k in r))
    lines = ["\t".join(cols)] + ["\t".join(str(r.get(c, "") if r.get(c) is not None else "").replace("\t", " ")
                                           for c in cols) for r in rows]
    return Response("﻿" + "\n".join(lines), media_type="text/tab-separated-values; charset=utf-8",
                    headers={"Content-Disposition": f"attachment; filename={filename}"})


def start_job(fn, label: str = "job", kind: str = "job") -> JSONResponse:
    """Run fn(progress) on the server as a tracked job (jobs.py); the page gets the job id and polls /api/job/<id>."""
    return JSONResponse({"job": jobs.start(fn, label, kind)})


CS_COLS = [("kind", "종류", True, True), ("decision", "결정", True, True), ("integrity", "정합성", True, True), ("seq", "서열 검증", True, True),
           ("conflict_fields", "불일치 필드", True, False), ("chosen_fields", "선택한 필드", False, False), ("accession", "Accession", True, False),
           ("accession_version", "Ver", True, False), ("genome_name", "Genome Name", True, False),
           ("completeness", "Completeness", True, True), ("genome_status", "NCBI Status", False, True), ("genome_length", "Length", True, False),
           ("strain", "Strain", True, False), ("host_name", "Host", True, False),
           ("clade", "Clade", None, False), ("lineage", "Lineage", None, False), ("serotype", "Serotype", None, False),
           ("genotype", "Genotype", None, False), ("subtype", "Subtype", None, False), ("outbreak", "Outbreak", None, False),
           ("isolation_country", "Country", True, False), ("geographic_location", "Geographic Location", False, False),
           ("collection_date", "Collection Date", True, False), ("isolation_source", "Isolation Source", False, False),
           ("sra_accession", "SRA Accession", True, False), ("bioproject_accession", "BioProject", False, False),
           ("biosample_accession", "BioSample", False, False), ("sequencing_platform", "Sequencing Platform", False, False),
           ("assembly_method", "Assembly Method", False, False), ("has_sra", "SRA", False, False),
           ("submitter", "Submitter", False, False), ("bvbrc_genome_id", "BV-BRC ID", False, False),
           ("release_date", "Release Date", False, False)]
CS_FACETS = [("kind", "종류"), ("decision", "결정"), ("integrity", "정합성(소스 간 일치)"), ("seq", "서열 검증"), ("completeness", "Completeness"), ("host_name", "Host"),
             ("isolation_country", "Country"), ("collection_year", "Collection Year"), ("has_sra", "SRA")]
GN_COLS = [("genome_name", "Genome Name", True), ("strain", "Strain", True), ("genbank", "GenBank Accessions", True),
           ("genome_length", "Size", True), ("completeness", "Completeness", True, True), ("ref_coverage", "Ref Coverage %", False),
           ("genome_status", "NCBI Status", False, True), ("collection_date", "Collection Date", True),
           ("isolation_country", "Isolation Country", True), ("host_name", "Host Name", True),
           ("clade", "Clade", None), ("lineage", "Lineage", None), ("serotype", "Serotype", None),
           ("genotype", "Genotype", None), ("subtype", "Subtype", None), ("outbreak", "Outbreak", None),
           ("geographic_location", "Geographic Location", False), ("isolation_source", "Isolation Source", False),
           ("sra_accession", "SRA Accession", True), ("bioproject_accession", "BioProject", False),
           ("biosample_accession", "BioSample", False), ("sequencing_platform", "Sequencing Platform", False),
           ("assembly_method", "Assembly Method", False), ("has_sra", "SRA Available", False), ("has_fasta", "FASTA", False),
           ("genome_quality", "Genome Quality", False), ("genome_id", "Genome ID", False), ("bvbrc_genome_id", "BV-BRC ID", False),
           ("qc", "QC", True, True), ("nc_clade", "Nextclade Clade", None), ("nc_lineage", "Nextclade Lineage", None),
           ("gc_content", "GC %", False), ("release_date", "Release Date", False), ("date_modified", "Modified", False)]
TAX_FACETS = [("clade", "Clade"), ("lineage", "Lineage"), ("serotype", "Serotype"), ("genotype", "Genotype"), ("subtype", "Subtype")]
GN_FACETS = [("completeness", "Completeness"), ("genome_status", "NCBI Status"), ("host_name", "Host"), ("isolation_country", "Country"),
             ("collection_year", "Collection Year"), ("has_sra", "SRA"), ("has_fasta", "FASTA"), ("qc", "QC"), ("genome_quality", "Quality")]


def _cols(cols, rows):
    """Optional taxonomy columns (on=None) are shown only when at least one row has a value."""
    out = []
    for c in cols:
        k, label, on = c[0], c[1], c[2]
        d = {"k": k, "label": label, "on": any(r.get(k) for r in rows) if on is None else on}
        if len(c) > 3:
            d["cls"] = c[3]
        out.append(d)
    return out


def _facets(facets, rows):
    extra = [f for f in TAX_FACETS if any(r.get(f[0]) for r in rows)]
    base = list(facets)
    at = next((i for i, f in enumerate(base) if f[0] == "collection_year"), len(base))
    return [{"k": k, "label": l} for k, l in base[:at] + extra + base[at:]]


def pending_total() -> int:
    return _rows("SELECT COUNT(*) n FROM changeset WHERE status IN ('READY','BLOCKED')")[0]["n"]


def dashboard_data(g: list[dict], prov: list[dict], last_checked) -> dict:
    from collections import Counter

    def top(key, n=15):
        c = Counter((r.get(key) or "(없음)") for r in g)
        items = c.most_common(n)
        rest = sum(c.values()) - sum(v for _, v in items)
        return {"labels": [k for k, _ in items] + (["기타"] if rest else []),
                "values": [v for _, v in items] + ([rest] if rest else [])}

    def by_prefix(key, n, empty="(없음)"):
        c = Counter((r.get(key) or "")[:n] or empty for r in g)
        keys = sorted(k for k in c if k != empty)
        return {"labels": keys + ([empty] if empty in c else []),
                "values": [c[k] for k in keys] + ([c[empty]] if empty in c else [])}

    lens = [r["genome_length"] for r in g if r.get("genome_length")]
    hist = {"labels": [], "values": []}
    if lens:
        lo, hi = min(lens), max(lens)
        step = max(1, -(-(hi - lo + 1) // 20))
        bins = Counter((x - lo) // step for x in lens)
        for i in range(max(bins) + 1):
            hist["labels"].append(f"{lo + i * step:,}")
            hist["values"].append(bins.get(i, 0))
    fields = ["host_name", "isolation_country", "collection_date", "isolation_source", "strain",
              "geographic_location", "genotype", "subtype", "genome_quality", "publication",
              "sra_accession", "bioproject_accession", "biosample_accession", "clade", "lineage"]
    n = len(g) or 1
    fill = sorted(({"f": f, "pct": round(100 * sum(1 for r in g if r.get(f)) / n, 1)} for f in fields),
                  key=lambda x: x["pct"])
    status = Counter(r.get("completeness") or r.get("genome_status") for r in g)
    return {
        "total": len(g), "complete_pct": round(100 * status.get("Complete", 0) / n, 1),
        "near_pct": round(100 * status.get("Near complete", 0) / n, 1),
        "avg_len": round(sum(lens) / len(lens)) if lens else 0, "last_checked": last_checked,
        "pending": pending_total(), "prov": {p["source"]: p["n"] for p in prov},
        "monthly": by_prefix("release_date", 7), "year": by_prefix("collection_date", 4),
        "country": top("isolation_country"), "host": top("host_name", 8), "status": top("completeness", 6),
        "clade": top("clade", 10) if any(r.get("clade") for r in g) else None,
        "lineage": top("lineage", 15) if any(r.get("lineage") for r in g) else None,
        "serotype": top("serotype", 10) if any(r.get("serotype") for r in g) else None,
        "hist": hist, "fill": {"labels": [x["f"] for x in fill], "values": [x["pct"] for x in fill]},
    }


def _ppx_slug(watch_id: int, rows: list[dict]) -> str | None:
    """Pathoplexus organism of a watch; discovered once from a joined accession and then stored."""
    cur = _rows("SELECT pathoplexus_organism FROM watchlist WHERE id=?", watch_id)
    if cur and cur[0]["pathoplexus_organism"]:
        return cur[0]["pathoplexus_organism"]
    acc = next((r["pathoplexus_accession"] for r in rows if r.get("pathoplexus_accession")), None)
    if not acc:
        return None
    from . import sources
    slug = sources.pathoplexus_slug_for(acc)
    if slug:
        con = connect()
        try:
            con.execute("UPDATE watchlist SET pathoplexus_organism=? WHERE id=?", (slug, watch_id))
            con.commit()
        finally:
            con.close()
    return slug


def _item_ids(form) -> list[int]:
    """Selected changeset item ids. Sent as ONE comma-separated field (`items`): one field per id would hit
    Starlette's 1000-field form limit as soon as more than 1000 rows are selected. Legacy `item` fields still work."""
    ids = [int(x) for x in form.getlist("item") if str(x).isdigit()]
    ids += [int(x) for x in str(form.get("items", "")).split(",") if x.strip().isdigit()]
    return ids


def create_app(scheduler: bool = True) -> FastAPI:
    log.setup()
    jobs.mark_interrupted()      # whatever the table still calls running died with the previous process
    app = FastAPI(title="Virus Genome Metadata Curator")

    def pending_count():
        return pending_total()

    def current_watch(request: Request, watch: int, watches: list[dict]) -> int:
        """Explicit ?watch= wins, then the cookie, then the first virus that has genomes."""
        ids = {x["id"] for x in watches}
        if watch in ids:
            return watch
        try:
            c = int(request.cookies.get("gc_watch", "0"))
        except ValueError:
            c = 0
        if c in ids:
            return c
        return next((x["id"] for x in watches if x["n"]), watches[0]["id"] if watches else 0)

    def schema_alert():
        """Latest schema check says something needs attention? Cached for a minute (it is read on every page)."""
        if time.time() - _alert_cache["t"] > 60:
            try:
                from . import schema_check
                a = schema_check.alert()
                _alert_cache.update(t=time.time(), v={"verdict": a["verdict"], "at": a["checked_at"]} if a else None, tl=schema_check.traffic())
            except Exception:
                _log.warning("could not read the schema alert", exc_info=True)
                _alert_cache.update(t=time.time(), v=None, tl=None)
        return _alert_cache["v"]

    def render(request, name, **ctx):
        resp = T.TemplateResponse(request, name, {"badge": pending_count(), "schema_alert": schema_alert(), "traffic": _alert_cache.get("tl"), "app_version": __version__, "watch": ctx.get("watch") or request.cookies.get("gc_watch", ""), **ctx})
        if ctx.get("watch"):
            resp.set_cookie("gc_watch", str(ctx["watch"]), max_age=86400 * 365, samesite="lax")
        return resp

    @app.get("/")
    def index(request: Request, msg: str = ""):
        return render(request, "index.html", msg=msg, nav="cur",
                      watches=_rows("SELECT w.*, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n_genomes FROM watchlist w"),
                      changesets=_changeset_cards())

    STATUS_TEXT = {"PENDING": ("준비 중", "warn"), "READY": ("검토 대기", "warn"), "BLOCKED": ("적용 불가", "bad"), "APPLIED": ("적용 완료", "ok")}

    def _changeset_cards():
        out = []
        for c in _rows("SELECT c.*, w.name wname FROM changeset c JOIN watchlist w ON w.id=c.watch_id "
                       "WHERE c.status!='DISCARDED' ORDER BY c.id DESC LIMIT 30"):
            sm = json.loads(c["summary_json"] or "{}")
            dec = {r["decision"]: r["n"] for r in _rows("SELECT decision, COUNT(*) n FROM changeset_item WHERE changeset_id=? GROUP BY decision", c["id"])}
            n_conf = 0
            if c["status"] != "APPLIED":
                ag = curate.source_agreement(c["watch_id"])
                n_conf = sum(1 for r in _rows("SELECT accession FROM changeset_item WHERE changeset_id=?", c["id"])
                             if ag.get(r["accession"], {}).get("conflict"))
            label, tone = STATUS_TEXT.get(c["status"], (c["status"], "warn"))
            pend, appr = dec.get("pending", 0), dec.get("approved", 0)
            if c["status"] == "APPLIED":
                nxt = "모두 반영되었습니다."
            elif c["status"] == "BLOCKED":
                nxt = "Dry-run에서 문제가 발견되어 적용할 수 없습니다. 검토 화면에서 원인을 확인하세요."
            elif pend:
                nxt = f"아직 결정하지 않은 항목이 {pend:,}건 있습니다. 검토하고 승인하세요."
            elif appr:
                nxt = f"승인된 {appr:,}건을 적용하면 내 데이터에 반영됩니다."
            else:
                nxt = "남은 작업이 없습니다."
            out.append({"id": c["id"], "wname": c["wname"], "created": (c["created"] or "")[:16], "status": c["status"], "label": label,
                        "tone": tone, "new": sm.get("new", 0), "changed": sm.get("changed", 0), "removed": sm.get("removed", 0),
                        "pending": pend, "approved": appr, "rejected": dec.get("rejected", 0), "hold": dec.get("hold", 0),
                        "applied": dec.get("applied", 0), "total": sum(dec.values()), "n_conf": n_conf, "next": nxt})
        return out

    def _same_origin(request: Request) -> bool:
        from urllib.parse import urlparse
        origin = request.headers.get("origin")
        return not origin or urlparse(origin).netloc == request.headers.get("host")

    @app.get("/api/watch/{wid}/impact")
    def watch_impact_api(wid: int):
        try:
            return curate.watch_impact(wid)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 404)

    @app.post("/watch/{wid}/update")
    def watch_update(request: Request, wid: int, name: str = Form(""), taxon: str = Form(""), enabled: str = Form("")):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        try:
            notes = curate.update_watch(wid, name, taxon, enabled == "1")
            msg = "저장했습니다: " + "; ".join(notes) if notes else "변경 사항이 없습니다."
        except Exception as e:
            _log.warning("watch update failed: %s: %s", type(e).__name__, e)
            msg = f"수정 실패: {e}"
        return RedirectResponse(f"/?msg={msg}", 303)

    @app.post("/watch/{wid}/delete")
    def watch_delete(request: Request, wid: int, confirm_name: str = Form("")):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        try:
            w = curate.watch_impact(wid)
            if confirm_name.strip() != w["name"]:
                msg = "삭제하지 않았습니다: 확인용 이름이 일치하지 않습니다."
            else:
                r = curate.delete_watch(wid)
                msg = (f"'{r['name']}' 삭제 완료 — 유전체 {r['genomes']}건, FASTA 파일 {r['files_removed']}개, "
                       f"changeset {r['changesets']}개를 지웠습니다.")
        except ValueError as e:
            msg = f"삭제 실패: {e}"
        return RedirectResponse(f"/?msg={msg}", 303)

    @app.post("/api/lookup/start")
    def api_lookup_start(request: Request, q: str = ""):
        """The same lookup as a server job, so it is not lost when the person clicks another menu; the result is stored with the job."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from .lookup import lookup_virus
        query = q.strip()
        if not query:
            return JSONResponse({"error": "이름을 입력하세요."}, 400)

        def work(p):
            p(5, "ICTV → NCBI Virus → BV-BRC → Pathoplexus 순서로 확인 중…")
            return {"redirect": None, "result": lookup_virus(query)}
        return start_job(work, f"바이러스 조회 · {query}", "lookup")

    @app.get("/api/lookup")
    def api_lookup(q: str = ""):
        """Check a typed virus name against ICTV, NCBI Virus, BV-BRC and Pathoplexus (read-only)."""
        from .lookup import lookup_virus
        try:
            return lookup_virus(q)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.post("/watch")
    def watch_add(request: Request, name: str = Form(...), taxon: str = Form("")):
        """Add a virus as a server job: it keeps going (and its outcome stays on the Jobs tab) when the person leaves the page."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        nm, tx = name.strip(), taxon.strip() or None

        def work(p):
            p(10, f"'{nm}'을(를) NCBI에서 바이러스로 확인하는 중…")
            curate.add_watch(nm, tx)
            p(95, "워치리스트에 등록했습니다")
            return f"/?msg=추가했습니다: {nm}"
        return start_job(work, f"바이러스 추가 · {nm}", "watch_add")

    @app.post("/watch/{name}/check")
    def check_now(name: str):
        def work(p):
            cs = curate.check(name, progress=p)
            if not cs:
                return "/?msg=변경 사항 없음"
            curate.dryrun(cs, progress=p)
            return f"/changeset/{cs}"
        return start_job(work, f"확인 · {name}", "check")

    @app.get("/api/version")
    def version_api():
        from .db import version_info
        return version_info()

    @app.get("/api/changeset/{cid}/item/{iid}/resolve")
    def resolve_view(cid: int, iid: int):
        from . import resolve
        con = connect()
        try:
            return resolve.item_view(con, cid, iid)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 404)
        finally:
            con.close()

    @app.post("/api/changeset/{cid}/item/{iid}/resolve")
    async def resolve_save(request: Request, cid: int, iid: int):
        """Save a person's choices between the sources' values for one item (choice source null = back to the default)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from . import resolve
        body = await request.json()
        con = connect()
        try:
            it = con.execute("SELECT i.accession, c.watch_id FROM changeset_item i JOIN changeset c ON c.id=i.changeset_id WHERE i.id=? AND c.id=?", (iid, cid)).fetchone()
            if not it:
                return JSONResponse({"error": "항목을 찾을 수 없습니다."}, 404)
            con.execute("BEGIN")
            res = resolve.save(con, it["watch_id"], it["accession"], body.get("choices", []))
            con.commit()
            _log.info("resolve %s: %s", it["accession"], res)
            return {**res, "view": resolve.item_view(con, cid, iid)}
        except ValueError as e:
            con.rollback()
            return JSONResponse({"error": str(e)}, 400)
        finally:
            con.close()

    @app.post("/api/changeset/{cid}/resolve_auto/preview")
    async def resolve_auto_preview_selected(request: Request, cid: int):
        """Preview for the selected rows only: body {"items": [item ids], "minimum": "medium"|"low"}."""
        from . import resolve
        body = await request.json()
        con = connect()
        try:
            m = body.get("minimum")
            pv = resolve.preview_auto(con, cid, m if m in resolve.CONF_RANK else "medium", [int(i) for i in body.get("items", [])])
            pv.pop("per_acc", None)
            return pv
        except (ValueError, TypeError) as e:
            return JSONResponse({"error": str(e)}, 404)
        finally:
            con.close()

    @app.get("/api/changeset/{cid}/resolve_auto")
    def resolve_auto_preview(cid: int, minimum: str = "medium"):
        from . import resolve
        con = connect()
        try:
            pv = resolve.preview_auto(con, cid, minimum if minimum in resolve.CONF_RANK else "medium")
            pv.pop("per_acc", None)
            return pv
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 404)
        finally:
            con.close()

    @app.post("/changeset/{cid}/resolve_auto")
    async def resolve_auto_apply(request: Request, cid: int, minimum: str = "medium"):
        """Choose the suggested value for every open disagreement of at least the given confidence."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from . import resolve
        try:
            body = await request.json()                  # optional: {"items": [item ids], "minimum": ...} for the selected rows only
        except Exception:
            body = {}
        if body.get("minimum") in resolve.CONF_RANK:
            minimum = body["minimum"]
        items = [int(i) for i in body["items"]] if isinstance(body.get("items"), list) else None
        con = connect()
        try:
            # no explicit BEGIN: the preview opens another connection (seqcheck.verdicts) and that must happen before this one writes,
            # otherwise SQLite refuses the write ("database is locked": the read snapshot of an open transaction is already stale)
            res = resolve.apply_auto(con, cid, minimum if minimum in resolve.CONF_RANK else "medium", items)
            con.commit()
            _log.info("resolve_auto changeset %s: %s", cid, res)
            return res
        except ValueError as e:
            con.rollback()
            return JSONResponse({"error": str(e)}, 400)
        finally:
            con.close()

    @app.get("/api/schema/traffic")
    def api_schema_traffic():
        from . import schema_check
        return schema_check.traffic()

    @app.get("/api/job/{jid}")
    def job(jid: str):
        return jobs.get(jid) or JSONResponse({"error": "unknown job", "done": True}, 404)

    @app.get("/api/jobs")
    def api_jobs(limit: int = 100, status: str = "", kind: str = ""):
        return jobs.list_jobs(max(1, min(limit, 300)), status, kind)

    @app.get("/api/jobs/summary")
    def api_jobs_summary(since: str = "", fail_since: str = ""):
        return jobs.summary(since, fail_since)

    def _sched_info(request: Request):
        """The scheduler's jobs: when each runs next and how the last run ended (None when the scheduler is off)."""
        s = getattr(request.app.state, "sched", None)
        if s is None:
            return None
        names = {"check": ("자동 확인", f"등록된 모든 바이러스를 {config.CHECK_INTERVAL_HOURS:g}시간마다 확인"),
                 "schema_daily": ("매일 스키마 점검", f"매일 {config.SCHEMA_CHECK_HOUR:02d}시, 바이러스마다 소스당 {config.SCHEMA_CHECK_SAMPLE}건")}
        out = []
        for j in s.get_jobs():
            title, what = names.get(j.id, (j.id, ""))
            last = _rows("SELECT status, started, finished, error FROM job WHERE kind='scheduled' AND label LIKE ? ORDER BY started DESC LIMIT 1", title + "%")
            out.append({"title": title, "what": what, "next": j.next_run_time.strftime("%Y-%m-%d %H:%M:%S") if j.next_run_time else None,
                        "last": last[0] if last else None})
        return out

    # ------------------------------------------------------------------------------------------------ export for other programs
    def _xf(request: Request) -> dict:
        qp = request.query_params
        lst = lambda k: [x.strip() for x in qp.get(k, "").split(",") if x.strip()] or None
        return {"require_sequence": qp.get("require_sequence", "1") not in ("0", "false"), "min_length": int(qp["min_length"]) if qp.get("min_length", "").isdigit() else None,
                "completeness_in": lst("completeness"), "qc_not": lst("qc_not"), "clade": qp.get("clade") or None}

    def _export_watch(watch: int):
        w = _rows("SELECT id, name FROM watchlist WHERE id=?", watch)
        return w[0] if w else None

    @app.get("/export")
    def export_page(request: Request, watch: int = 0):
        watches = _rows("SELECT w.id, w.name, (SELECT COUNT(*) FROM genome_sequence s JOIN genome g ON g.genome_id=s.genome_id WHERE g.watch_id=w.id) ns FROM watchlist w ORDER BY w.name")
        return render(request, "export.html", nav="exp", watches=watches, profiles={k: v["label"] for k, v in export.PROFILES.items()},
                      exports=export.list_exports(), exports_dir=str(export.exports_dir()))

    @app.get("/api/export/profiles")
    def api_export_profiles():
        return {k: {"label": v["label"], "id_column": v["id"], "columns": v["columns"]} for k, v in export.PROFILES.items()}

    @app.get("/api/export/{watch}/preview")
    def api_export_preview(request: Request, watch: int, profile: str = "nextstrain"):
        if not _export_watch(watch):
            return JSONResponse({"error": "unknown virus"}, 404)
        try:
            return export.preview(watch, profile, **_xf(request))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/export/{watch}/metadata.tsv")
    def api_export_metadata(request: Request, watch: int, profile: str = "nextstrain"):
        """The metadata table as TSV - the same rows, in the same order, as sequences.fasta."""
        w = _export_watch(watch)
        if not w:
            return JSONResponse({"error": "unknown virus"}, 404)
        try:
            cols, rows, _, _ = export.table(watch, profile, **_xf(request))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        return StreamingResponse(export.iter_metadata(cols, rows), media_type="text/tab-separated-values; charset=utf-8",
                                 headers={"Content-Disposition": f'attachment; filename="{export.slug(w["name"])}_metadata.tsv"'})

    @app.get("/api/export/{watch}/sequences.fasta")
    def api_export_fasta(request: Request, watch: int, profile: str = "nextstrain"):
        """FASTA records whose ids equal the profile's id column in metadata.tsv."""
        w = _export_watch(watch)
        if not w:
            return JSONResponse({"error": "unknown virus"}, 404)
        try:
            _, _, sel, _ = export.table(watch, profile, **_xf(request))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        return StreamingResponse(export.iter_fasta(profile, sel), media_type="text/plain; charset=utf-8",
                                 headers={"Content-Disposition": f'attachment; filename="{export.slug(w["name"])}_sequences.fasta"'})

    @app.post("/api/export/{watch}/build")
    def api_export_build(request: Request, watch: int, profile: str = "nextstrain"):
        """Write metadata.tsv, sequences.fasta, manifest.json and nextstrain-config.yaml into data/exports/<virus>/<timestamp>/ (a server job)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        w = _export_watch(watch)
        if not w:
            return JSONResponse({"error": "unknown virus"}, 404)
        if profile not in export.PROFILES:
            return JSONResponse({"error": "unknown profile"}, 400)
        flt = _xf(request)
        return start_job(lambda p: {"redirect": "/export", "result": export.write(watch, None, profile, p, **flt)}, f"내보내기 · {w['name']} ({profile})", "export")

    @app.get("/api/exports")
    def api_exports(virus: str = ""):
        return export.list_exports(virus or None)

    # ------------------------------------------------------------------------------------------------ QC (sequence stats, Nextclade)
    @app.get("/qc")
    def qc_page(request: Request, watch: int = 0):
        watches = _rows("SELECT w.id, w.name, w.taxon, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n, "
                       "(SELECT COUNT(*) FROM genome_sequence s JOIN genome g ON g.genome_id=s.genome_id WHERE g.watch_id=w.id) ns FROM watchlist w ORDER BY w.name")
        return render(request, "qc.html", nav="qc", watches=watches, tool_version=qc.nextclade_version(), tool_ok=qc.nextclade_bin() is not None,
                      summary={w["id"]: qc.summary(w["id"]) for w in watches},
                      datasets_map={w["id"]: qc.dataset_for(w["taxon"]) for w in watches})

    @app.get("/api/qc/datasets")
    def api_qc_datasets():
        return {"tool": qc.nextclade_version(), "rows": qc.dataset_status()}

    @app.post("/api/qc/run")
    def api_qc_run(request: Request, watch: int, force: int = 0, then_check: int = 0):
        """QC of one virus as a server job; then_check also runs 'check now' so Nextclade's clade/lineage enter the review screen."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        w = _rows("SELECT name FROM watchlist WHERE id=?", watch)
        if not w:
            return JSONResponse({"error": "unknown virus"}, 404)
        name = w[0]["name"]

        def work(p):
            res = qc.run(watch, lambda pct, m: p(pct * (0.6 if then_check else 1), m), bool(force))
            if then_check:
                p(62, "QC 결과를 소스와 함께 대조하는 중 (지금 확인)…")
                cs = curate.check(name, progress=lambda pct, m: p(62 + pct * 0.3, m))
                if cs:
                    curate.dryrun(cs, progress=lambda pct, m: p(92 + pct * 0.08, m))
                return {"redirect": f"/changeset/{cs}" if cs else "/qc", "result": res}
            return {"redirect": f"/qc?watch={watch}", "result": res}
        return start_job(work, f"QC · {name}", "qc")

    @app.post("/api/qc/dataset")
    def api_qc_dataset(request: Request, path: str):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        known = {qc.dataset_for(w["taxon"]) for w in _rows("SELECT taxon FROM watchlist")}
        if path not in known:
            return JSONResponse({"error": "등록된 바이러스에 연결된 데이터셋이 아닙니다."}, 400)
        return start_job(lambda p: {"redirect": "/qc", "result": qc.dataset_update(path, p)}, f"Nextclade 데이터셋 업데이트 · {path}", "qc")

    @app.get("/analysis")
    def analysis_page(request: Request):
        return render(request, "analysis.html", nav="ana", st=analysis.status())

    @app.get("/api/superset/status")
    def api_superset_status():
        return analysis.status()

    @app.post("/api/superset/start")
    def api_superset_start(request: Request):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        return start_job(lambda p: {"redirect": "/analysis", "result": {"url": analysis.start(p)}}, "Superset 켜기", "superset")

    @app.post("/api/superset/stop")
    def api_superset_stop(request: Request):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        return {"stopped": analysis.stop()}

    @app.get("/jobs")
    def jobs_page(request: Request):
        return render(request, "jobs.html", nav="job", sched=_sched_info(request), kinds=jobs.KINDS, scheduler_on=getattr(request.app.state, "sched", None) is not None)

    @app.get("/api/changeset/{cid}/impact")
    def changeset_impact_api(cid: int):
        try:
            return curate.changeset_impact(cid)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 404)

    @app.post("/changeset/{cid}/delete")
    def changeset_delete(request: Request, cid: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        try:
            r = curate.delete_changeset(cid)
            msg = f"changeset {cid}번 삭제 완료 — 항목 {r['items']}건 (이미 반영된 유전체는 그대로입니다)."
        except ValueError as e:
            msg = f"삭제 실패: {e}"
        return RedirectResponse(f"/?msg={msg}", 303)

    @app.get("/changeset/{cid}")
    def changeset(request: Request, cid: int):
        cs = _rows("SELECT c.*, w.name wname FROM changeset c JOIN watchlist w ON w.id=c.watch_id WHERE c.id=?", cid)[0]
        rows = []
        _c = connect()
        try:
            _ref = completeness.reference_map(_c).get(cs["watch_id"])
        finally:
            _c.close()
        agree = curate.source_agreement(cs["watch_id"])
        for i in _rows("SELECT * FROM changeset_item WHERE changeset_id=? ORDER BY kind, accession", cid):
            r = {k: v for k, v in json.loads(i["new_record_json"] or "{}").items()
                 if v is not None and k not in ("taxon_lineage_ids", "taxon_lineage_names")}
            if isinstance(r.get("submitter"), str):
                r["submitter"] = r["submitter"][:120]
            ag = agree.get(i["accession"], {})
            r["completeness"], r["ref_coverage"] = completeness.classify(r.get("genome_status"), r.get("genome_length"), _ref and _ref[1])
            r.update(id=i["id"], kind=i["kind"], decision=i["decision"], accession=i["accession"],
                     has_sra="Yes" if r.get("sra_accession") else "No",
                     integrity=("불일치" if ag.get("conflict") else "표기차이" if ag.get("cosmetic") else "선택 완료" if ag.get("chosen") else "일치"),
                     chosen_fields=", ".join(sorted(ag.get("chosen", []))),
                     seq={"same": "서열 동일", "diff": "서열 다름", "other": "확인 불가"}.get(ag.get("seq"), "미검증"),
                     conflict_fields=", ".join(sorted(ag.get("conflict", []))))
            if i["kind"] == "changed":
                r["_diff"] = json.loads(i["diff_json"] or "{}")
            rows.append(r)
        data = {"rows": rows, "selectable": True, "name": f"changeset_{cid}",
                "cols": _cols(CS_COLS, rows), "facets": _facets(CS_FACETS, rows)}
        n_conf = sum(1 for r in rows if r["integrity"] == "불일치")
        return render(request, "changeset.html", cs=cs, data=data, nav="cur", watch=cs["watch_id"], n_conf=n_conf,
                      summary=json.loads(cs["summary_json"] or "{}"),
                      dry=json.loads(cs["dryrun_json"] or "{}"))

    @app.post("/changeset/{cid}/decide")
    async def decide(request: Request, cid: int):
        form = await request.form()
        ids = _item_ids(form)
        action = form.get("decision")
        if action in ("approved", "rejected", "hold"):
            curate.decide(ids, action)
        return RedirectResponse(f"/changeset/{cid}", 303)

    def _apply_job(cid: int, fetch: bool = True):
        def work(p):
            res = curate.apply(cid, progress=p, fetch_sequences=fetch)
            msg = f"{res['applied']}건 적용 완료"
            if res.get("no_fasta"):
                msg += f" (신규 {res['no_fasta']}건은 메타데이터만 — 유전체 화면의 'Get FASTA'로 서열을 나중에 받을 수 있음)"
            if res.get("deferred"):
                msg += f", 서열이 바뀐 {res['deferred']}건은 서열이 필요해 승인 상태로 남김"
            if res.get("metadata_only"):
                msg += f" (메타데이터만 변경 {res['metadata_only']}건은 서열을 다시 받지 않음)"
            if res.get("failed"):
                first = next(iter(res["failed"].items()))
                msg += f", {len(res['failed'])}건 실패 (예: {first[0]} - {first[1]}). 실패 항목은 승인 상태로 남아 있으니 Apply를 다시 누르면 재시도합니다"
            return f"/history?msg={msg}"
        return start_job(work, f"적용 · changeset {cid}", "apply")

    @app.post("/changeset/{cid}/approve_apply")
    async def approve_apply(request: Request, cid: int):
        form = await request.form()
        curate.decide(_item_ids(form), "approved")
        return _apply_job(cid, form.get("fetch_sequences", "1") != "0")

    @app.post("/changeset/{cid}/apply")
    async def do_apply(request: Request, cid: int):
        form = await request.form()
        return _apply_job(cid, form.get("fetch_sequences", "1") != "0")

    @app.post("/genomes/sequences")
    async def genomes_sequences(request: Request):
        """Fetch FASTA for genomes that were applied metadata-only."""
        form = await request.form()
        watch = int(form.get("watch", 0))
        ids = [i for i in str(form.get("ids", "")).split(",") if i]

        def work(p):
            res = curate.fetch_missing_sequences(watch, ids or None, progress=p)
            msg = f"서열 {res['fetched']}건 다운로드 완료"
            if res["failed"]:
                first = next(iter(res["failed"].items()))
                msg += f", {len(res['failed'])}건 실패 (예: {first[0]} - {first[1]})"
            return f"/genomes?watch={watch}&msg={msg}"
        return start_job(work, "서열(FASTA) 받기", "fasta")

    @app.post("/changeset/{cid}/dryrun")
    def rerun(cid: int):
        curate.dryrun(cid)
        return RedirectResponse(f"/changeset/{cid}", 303)

    @app.get("/genomes")
    def genomes(request: Request, watch: int = 0, msg: str = ""):
        watches = _rows("SELECT w.id, w.name, w.taxon_id, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n FROM watchlist w")
        watch = current_watch(request, watch, watches)
        cur = next((w for w in watches if w["id"] == watch), None)
        rows = _rows("SELECT g.*, s.md5, s.gc_content FROM genome g LEFT JOIN genome_sequence s ON s.accession=g.accession "
                     "WHERE g.watch_id=? ORDER BY g.genome_id", watch)
        prov: dict[str, dict] = {}
        for p in _rows("SELECT p.genome_id, p.field, p.source FROM field_provenance p JOIN genome g ON g.genome_id=p.genome_id "
                       "WHERE g.watch_id=?", watch):
            prov.setdefault(p["genome_id"], {})[p["field"]] = p["source"]
        lineage = (rows[0].get("taxon_lineage_names") or "").split(",") if rows else []
        _c = connect()
        try:
            ref = completeness.reference_map(_c).get(watch)
        finally:
            _c.close()
        qcm = qc.per_genome(watch)
        for n, r in enumerate(rows):
            q = qcm.get(r["accession"], {})
            r["qc"], r["nc_clade"], r["nc_lineage"] = q.get("qc", "미실행"), q.get("nc_clade"), q.get("nc_lineage")
            r.pop("taxon_lineage_ids", None)
            r.pop("taxon_lineage_names", None)
            r["id"] = n
            r["genbank"] = f"{r['accession']}.{r['accession_version']}" if r.get("accession_version") else r["accession"]
            r["_prov"] = prov.get(r["genome_id"], {})
            r["has_sra"] = "Yes" if r.get("sra_accession") else "No"
            r["has_fasta"] = "Yes" if r.get("md5") else "No"
            r["completeness"], r["ref_coverage"] = completeness.classify(r.get("genome_status"), r.get("genome_length"), ref and ref[1])
        ppx = _ppx_slug(watch, rows)
        data = {"rows": rows, "name": "genomes", "watch": cur and cur["name"], "lineage": lineage, "ppx": ppx,
                "ref_note": completeness.describe(ref),
                "cols": _cols(GN_COLS, rows), "facets": _facets(GN_FACETS, rows)}
        return render(request, "genomes.html", data=data, nav="gen", watches=watches, watch=watch, msg=msg)

    @app.post("/genomes/fasta")
    async def genomes_fasta(request: Request):
        form = await request.form()
        ids = [i for i in str(form.get("ids", "")).split(",") if i]
        out = []
        for r in _rows("SELECT file_path FROM genome_sequence WHERE genome_id IN (%s)" % ",".join("?" * len(ids)), *ids) if ids else []:
            try:
                out.append(config.data_path(r["file_path"]).read_text(encoding="utf-8"))
            except OSError:
                pass
        return Response("".join(out), media_type="text/plain; charset=utf-8",
                        headers={"Content-Disposition": "attachment; filename=genomes.fna"})

    @app.get("/genomes.tsv")
    def genomes_tsv():
        return _tsv(_rows("SELECT g.*, s.gc_content, s.md5, s.file_path FROM genome g "
                          "LEFT JOIN genome_sequence s ON s.accession=g.accession ORDER BY g.genome_id"), "genomes.tsv")

    @app.post("/integrity/seqplan")
    def integrity_seqplan(request: Request, watch: int = Form(0)):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from . import seqcheck
        return start_job(lambda p: (seqcheck.build_plan(watch, None, progress=p), f"/integrity?watch={watch}#seq")[1], "서열 검증 계획", "seqcheck")

    @app.post("/integrity/seqcheck")
    def integrity_seqcheck(request: Request, watch: int = Form(0)):
        """Execute the saved plan (the user has seen and confirmed it)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from . import seqcheck
        return start_job(lambda p: (seqcheck.execute(watch, progress=p), f"/integrity?watch={watch}#seq")[1], "서열 검증 실행", "seqcheck")

    @app.post("/changeset/{cid}/seqcheck")
    async def changeset_seqcheck(request: Request, cid: int):
        """Plan (dry-run) for the selected items, then send the user to the plan page to confirm."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from . import seqcheck
        form = await request.form()
        ids = _item_ids(form)
        cs = _rows("SELECT watch_id FROM changeset WHERE id=?", cid)
        if not cs or not ids:
            return JSONResponse({"error": "선택된 항목이 없습니다."}, 400)
        accs = [r["accession"] for r in _rows("SELECT accession FROM changeset_item WHERE changeset_id=? AND id IN (%s)" % ",".join("?" * len(ids)), cid, *ids)]
        wid = cs[0]["watch_id"]
        return start_job(lambda p: (seqcheck.build_plan(wid, accs, progress=p), f"/integrity?watch={wid}#seq")[1], "선택 항목 서열 검증 계획", "seqcheck")

    @app.get("/schema")
    def schema_page(request: Request, msg: str = ""):
        from . import schema_check
        return render(request, "schema.html", nav="sch", msg=msg, rep=schema_check.latest(), daily=schema_check.latest_daily(),
                      label=schema_check.LABEL, sample=config.SCHEMA_CHECK_SAMPLE, hour=config.SCHEMA_CHECK_HOUR)

    @app.post("/schema/check")
    def schema_run(request: Request):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from . import schema_check
        return start_job(lambda p: (schema_check.run(None, progress=p), "/schema")[1], "스키마 점검", "schema")

    @app.post("/schema/daily")
    def schema_daily(request: Request):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from . import schema_check

        def work(p):
            r = schema_check.run_light(progress=p)
            _alert_cache["t"] = 0
            return "/schema"
        return start_job(work, "바이러스별 스키마 점검", "schema")

    @app.post("/schema/accept")
    def schema_accept(request: Request):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        from . import schema_check
        try:
            msg = "새 기준으로 저장했습니다: " + ", ".join(schema_check.LABEL.get(x, x) for x in schema_check.accept())
            _alert_cache["t"] = 0
        except ValueError as e:
            msg = str(e)
        return RedirectResponse(f"/schema?msg={msg}", 303)

    @app.get("/integrity")
    def integrity(request: Request, watch: int = 0, field: str = "", kind: str = "conflict", files: int = 0):
        from . import reconcile, seqcheck, verify
        watches = _rows("SELECT id, name FROM watchlist ORDER BY name")
        wid = watch or (int(request.cookies.get("gc_watch") or 0) if request.cookies.get("gc_watch", "").isdigit() else 0)
        if wid not in {w["id"] for w in watches}:
            wid = watches[0]["id"] if watches else 0
        rep = reconcile.load(wid) if wid else None
        tot = cmp_ = 0
        for byp in (rep or {}).get("fields", {}).values():
            for d in byp.values():
                tot += d["match"] + d["cosmetic"] + d["conflict"]
                cmp_ += d["match"] + d["cosmetic"]
        pairs = [{"key": f"{a}|{b}", "label": f"{reconcile.LABEL[a]} ↔ {reconcile.LABEL[b]}"}
                 for a, b in (("ncbi", "bvbrc"), ("ncbi", "pathoplexus"), ("bvbrc", "pathoplexus"))]
        return render(request, "integrity.html", nav="int", watch=wid, wid=wid, watches=watches, rep=rep, pairs=pairs,
                      label=reconcile.LABEL, agree=100 * cmp_ / tot if tot else 100.0, tot_cmp=tot,
                      ext=sum(rep["coverage"]["not_in_ncbi"].values()) if rep else 0,
                      field=field, kind=kind, examples=reconcile.conflicts(wid, field, kind) if field and wid else [],
                      checks=verify.run(wid or None, check_files=bool(files)), files=files,
                      seq=seqcheck.summary(wid) if wid else None, plan=seqcheck.load_plan(wid) if wid else None, seq_rows=seqcheck.results(wid) if wid else [])

    @app.get("/dashboard")
    def dashboard(request: Request, watch: int = 0):
        watches = _rows("SELECT w.id, w.name, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n FROM watchlist w")
        watch = current_watch(request, watch, watches)
        g = _rows("SELECT * FROM genome WHERE watch_id=?", watch)
        _c = connect()
        try:
            _ref = completeness.reference_map(_c).get(watch)
        finally:
            _c.close()
        for r in g:
            r["completeness"], _ = completeness.classify(r.get("genome_status"), r.get("genome_length"), _ref and _ref[1])
        prov = _rows("SELECT p.source, COUNT(*) n FROM field_provenance p JOIN genome g ON g.genome_id=p.genome_id "
                     "WHERE g.watch_id=? GROUP BY p.source", watch)
        last = _rows("SELECT last_checked FROM watchlist WHERE id=?", watch)
        return render(request, "dashboard.html", nav="dash", watches=watches, watch=watch,
                      data=dashboard_data(g, prov, last[0]["last_checked"] if last else None))

    @app.get("/api/deploy/preview")
    def deploy_preview():
        from . import deploy
        try:
            return deploy.preview()
        except Exception as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.post("/deploy")
    def deploy_now(request: Request):
        """Build docs/index.html from the local DB, commit it and push (GitHub Pages)."""
        from urllib.parse import urlparse

        from . import deploy
        origin = request.headers.get("origin")
        if origin and urlparse(origin).netloc != request.headers.get("host"):
            return JSONResponse({"error": "cross-site request refused"}, 403)

        def work(p):
            return "/?msg=" + deploy.publish(p)
        return start_job(work, "공개 페이지 배포", "deploy")

    @app.get("/public")
    def public_view():
        """The same read-only page as docs/index.html, built live from the local database."""
        from .static_site import build_html
        try:
            return HTMLResponse(build_html("https://github.com/CHOICHI-HWAN/GenomeCurator"))
        except ValueError as e:
            return HTMLResponse(f"<p>{e}</p>", status_code=404)

    @app.get("/history")
    def history(request: Request, msg: str = "", lvl: str = "", n: int = 200, q: str = ""):
        n = max(20, min(n, 2000))
        return render(request, "history.html", msg=msg, nav="his", lvl=lvl, n=n, q=q, log_dir=str(log.log_dir()),
                      log=_rows("SELECT * FROM audit_log ORDER BY id DESC LIMIT 200"), srv=log.tail(n, lvl, q))

    if scheduler:
        from apscheduler.schedulers.background import BackgroundScheduler

        def sched_job():
            _log.info("scheduled check run (every %s h)", config.CHECK_INTERVAL_HOURS)
            for w in _rows("SELECT name FROM watchlist WHERE enabled=1"):
                def work(p, name=w["name"]):
                    cs = curate.check(name, progress=p)
                    if cs:
                        curate.dryrun(cs, progress=p)
                    return f"/changeset/{cs}" if cs else "/"
                jobs.run_sync(f"자동 확인 · {w['name']}", "scheduled", work)     # recorded like any job: status, log, error

        def schema_job():
            from . import schema_check
            _log.info("daily schema check start (%d records per source per virus)", config.SCHEMA_CHECK_SAMPLE)

            def work(p):
                rep = schema_check.run_light(progress=p)
                _alert_cache["t"] = 0   # show a new problem right away
                _log.info("daily schema check verdict: %s", rep["verdict"])
                return "/schema"
            jobs.run_sync("매일 스키마 점검", "scheduled", work)

        sched = BackgroundScheduler()
        sched.add_job(sched_job, "interval", hours=config.CHECK_INTERVAL_HOURS, id="check")
        sched.add_job(schema_job, "cron", hour=config.SCHEMA_CHECK_HOUR, minute=0, id="schema_daily", misfire_grace_time=3600, coalesce=True)
        sched.start()
        app.state.sched = sched
        for j in sched.get_jobs():
            _log.info("scheduled job %s: next run %s", j.id, j.next_run_time)
        app.router.on_shutdown.append(lambda: sched.shutdown(wait=False))
    return app
