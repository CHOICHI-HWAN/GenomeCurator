import json
import os
import time
from pathlib import Path

from fastapi import FastAPI, Form, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse, Response, StreamingResponse
from fastapi.templating import Jinja2Templates

from . import cluster, peartree, beastprep, cladedefs, clock, clockclass, continents, focus, hyphyexport, intake, lineagetree, mutgantt, nextcladecall, pipeline, support, variantcall, variantflow, __version__, analysis, completeness, datasets, exclusions, export, cooccur, geo, heatmap, vocab, insight, lineages, mutprofile, nextstrain, novel, outbreak, pathoplexus, private, qc, squirrel, config, curate, jobs, log
from .db import audit_now, connect

_log = log.get("web")
_alert_cache = {"t": 0.0, "v": None}
_STARTED = time.time()                    # when this server process started: code files newer than that are not in the running app yet
_stale_cache = {"t": 0.0, "v": False}


def code_is_newer_than_server() -> bool:
    """True when a source or template file changed after this process started (the app does not reload itself): the screen then tells the user to restart. Checked at most every 20 s."""
    if time.time() - _stale_cache["t"] > 20:
        try:
            root = Path(__file__).resolve().parent
            newest = max((f.stat().st_mtime for pat in ("*.py", "templates/*.html", "resources/vocab/*.tsv") for f in root.glob(pat)), default=0)
            _stale_cache.update(t=time.time(), v=newest > _STARTED + 1)
        except OSError:
            _stale_cache.update(t=time.time(), v=False)
    return _stale_cache["v"]

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
           ("biosample_accession", "BioSample", False, False), ("sequencing_platform", "Sequencing Platform", False, False), ("platform", "Platform", False, False),
           ("assembly_method", "Assembly Method", False, False), ("has_sra", "SRA", False, False),
           ("submitter", "Submitter", False, False), ("bvbrc_genome_id", "BV-BRC ID", False, False),
           ("release_date", "Release Date", False, False)]
CS_FACETS = [("kind", "종류"), ("decision", "결정"), ("integrity", "정합성(소스 간 일치)"), ("seq", "서열 검증"), ("completeness", "Completeness"), ("host_name", "Host"),
             ("isolation_country", "Country"), ("collection_year", "Collection Year"), ("has_sra", "SRA")]
GN_COLS = [("genome_name", "Genome Name", True), ("strain", "Strain", True), ("genbank", "GenBank Accessions", True),
           ("genome_length", "Size", True), ("completeness", "Completeness", True, True), ("excluded", "데이터셋", True, True), ("ref_coverage", "Ref Coverage %", False),
           ("genome_status", "NCBI Status", False, True), ("collection_date", "Collection Date", True),
           ("isolation_country", "Isolation Country", True), ("host_name", "Host Name", True),
           ("clade", "Clade", None), ("lineage", "Lineage", None), ("serotype", "Serotype", None),
           ("genotype", "Genotype", None), ("subtype", "Subtype", None), ("outbreak", "Outbreak", None),
           ("geographic_location", "Geographic Location", False), ("isolation_source", "Isolation Source", False),
           ("sra_accession", "SRA Accession", True), ("bioproject_accession", "BioProject", False),
           ("biosample_accession", "BioSample", False), ("sequencing_platform", "Sequencing Platform", False), ("platform", "Platform", False),
           ("assembly_method", "Assembly Method", False), ("has_sra", "SRA Available", False), ("has_fasta", "FASTA", False),
           ("genome_quality", "Genome Quality", False), ("genome_id", "Genome ID", False), ("bvbrc_genome_id", "BV-BRC ID", False),
           ("qc", "QC", True, True), ("nc_clade", "Nextclade Clade", None), ("nc_lineage", "Nextclade Lineage", None),
           ("gc_content", "GC %", False), ("release_date", "Release Date", False), ("date_modified", "Modified", False)]
TAX_FACETS = [("clade", "Clade"), ("lineage", "Lineage"), ("serotype", "Serotype"), ("genotype", "Genotype"), ("subtype", "Subtype")]
GN_FACETS = [("completeness", "Completeness"), ("excluded", "데이터셋 (제외 표시)"), ("genome_status", "NCBI Status"), ("host_name", "Host"), ("isolation_country", "Country"),
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
    nextstrain.mark_interrupted()
    app = FastAPI(title="Virus Genome Metadata Curator")

    @app.middleware("http")
    async def other_machines_reach_the_cluster_api_only(request: Request, call_next):
        """The app has no log-in: every page shows (and every route changes) the curated and private data. When it listens on the network for the cluster
        (`serve --host 0.0.0.0`), another machine may therefore use only the worker API, which asks for a worker key; everything else answers 403.
        GC_ALLOW_REMOTE_UI=1 lifts this (for a trusted network where the whole app is meant to be shared)."""
        host = request.client.host if request.client else ""
        if host not in ("127.0.0.1", "::1", "localhost", "testclient") and not request.url.path.startswith("/api/cluster/w/") and os.environ.get("GC_ALLOW_REMOTE_UI") != "1":
            return JSONResponse({"error": "this machine may use the cluster worker API only; open the app on the main server itself"}, 403)
        return await call_next(request)

    @app.middleware("http")
    async def refuse_cross_site_writes(request: Request, call_next):
        """The app runs on the user's own machine, so a web page opened in the same browser must not be able to make it change anything: a write (POST / PUT / PATCH / DELETE)
        whose Origin is another site is refused for every route, including the ones that predate the per-route check."""
        if request.method in ("POST", "PUT", "PATCH", "DELETE"):
            from urllib.parse import urlparse
            origin = request.headers.get("origin")
            if (origin and origin != "null" and urlparse(origin).netloc != request.headers.get("host")) or origin == "null" or request.headers.get("sec-fetch-site") == "cross-site":
                return JSONResponse({"error": "cross-site request refused"}, 403)
        return await call_next(request)

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
        resp = T.TemplateResponse(request, name, {"badge": pending_count(), "schema_alert": schema_alert(), "traffic": _alert_cache.get("tl"), "app_version": __version__, "stale_code": code_is_newer_than_server(), "watch": ctx.get("watch") or request.cookies.get("gc_watch", ""), **ctx})
        if ctx.get("watch"):
            resp.set_cookie("gc_watch", str(ctx["watch"]), max_age=86400 * 365, samesite="lax")
        return resp

    @app.get("/")
    def index(request: Request, msg: str = ""):
        ws = _rows("SELECT w.*, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id AND IFNULL(g.genome_status,'')!='Removed') n_genomes FROM watchlist w")
        for w in ws:
            w["support"] = support.for_virus(w["name"], w["taxon"])
        return render(request, "index.html", msg=msg, nav="cur", watches=ws,
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
            try:
                curate.add_watch(nm, tx)
            except curate.AlreadyWatched as e:                  # already there: say so and go to it, it is not an error
                return f"/genomes?watch={e.watch_id}&msg=이미 워치리스트에 있습니다: {e.name} (유전체 화면으로 이동했습니다)"
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

    def _ds(request: Request, watch: int):
        """?dataset=<id>: the named dataset decides who is in (its rules replace the other filters). Downloads never carry private genomes."""
        v = request.query_params.get("dataset", "")
        if not v:
            return None
        d = datasets.get(int(v)) if v.isdigit() else None
        if not d or d["watch_id"] != watch:
            raise ValueError("그 바이러스에 없는 데이터셋입니다.")
        return d

    def _export_watch(watch: int):
        w = _rows("SELECT id, name FROM watchlist WHERE id=?", watch)
        return w[0] if w else None

    @app.get("/export")
    def export_page(request: Request, watch: int = 0):
        watches = _rows("SELECT w.id, w.name, (SELECT COUNT(*) FROM genome_sequence s JOIN genome g ON g.genome_id=s.genome_id WHERE g.watch_id=w.id) ns FROM watchlist w ORDER BY w.name")
        return render(request, "export.html", nav="exp", watches=watches, profiles={k: v["label"] for k, v in export.PROFILES.items()},
                      exports=export.list_exports(), exports_dir=str(export.exports_dir()))

    @app.get("/api/exclusions/{watch}")
    def api_exclusions(watch: int, scope: str = ""):
        return exclusions.listing(watch, scope or None)

    @app.post("/api/exclusions/{watch}")
    async def api_exclusions_set(request: Request, watch: int):
        """Leave genomes out of the dataset (or put them back): body {ids: [...], scope: public|private, action: exclude|include, reason}. Nothing is deleted."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        body = await request.json()
        ids, scope = [str(i) for i in body.get("ids", [])], body.get("scope", "public")
        try:
            if body.get("action") == "include":
                return exclusions.remove(watch, ids, scope)
            return exclusions.add(watch, ids, scope, str(body.get("reason") or ""))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    # ------------------------------------------------------------------------------------------------ datasets
    @app.get("/datasets")
    def datasets_page(request: Request, watch: int = 0):
        watches = _rows("SELECT id, name FROM watchlist ORDER BY name")
        wid = watch or next((w["id"] for w in watches if "monkeypox" in w["name"].lower()), watches[0]["id"] if watches else 0)
        return render(request, "datasets.html", nav="ds", watches=watches, wid=wid, builds=datasets.BUILDS)

    @app.get("/api/datasets")
    def api_datasets(watch: int = 0):
        return datasets.listing(watch or None)

    @app.get("/api/datasets/facets")
    def api_datasets_facets(watch: int):
        try:
            return datasets.facets(watch)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 404)

    def _ds_from_body(body: dict) -> dict:
        return {"id": body.get("id"), "watch_id": int(body["watch"]), "name": body.get("name") or "", "description": body.get("description") or "",
                "rules": datasets.normalize(body.get("rules")), "include_private": bool(body.get("include_private")),
                "manual_include": [str(x).strip() for x in body.get("manual_include", []) if str(x).strip()],
                "manual_exclude": [str(x).strip() for x in body.get("manual_exclude", []) if str(x).strip()], "build": body.get("build") or ""}

    @app.post("/api/datasets/preview")
    async def api_datasets_preview(request: Request):
        """Counts for a dataset that is still being edited (nothing is saved): body = the dataset fields + watch."""
        body = await request.json()
        try:
            ds = _ds_from_body(body)
            pv = export.preview(ds["watch_id"], "nextstrain", include_private=ds["include_private"], dataset=ds)
            cols, rows, sel, _ = export.table(ds["watch_id"], "nextstrain", ds)
            from collections import Counter
            top = lambda k: Counter(r.get(k) or "(비어 있음)" for r in rows).most_common(8)
            pt, pm = len(export.private_table(ds["watch_id"])[1]), len(export.private_table(ds["watch_id"], ds)[1])
            return {**pv, "private_total": pt, "private_match": pm, "by_clade": top("clade"), "by_country": top("country"), "by_year": sorted(Counter((r["date"] or "????")[:4] for r in rows).items())[-12:],
                    "summary": datasets.describe(ds["rules"])}
        except (ValueError, KeyError, TypeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.post("/api/datasets")
    async def api_datasets_save(request: Request):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        body = await request.json()
        try:
            return datasets.save(int(body["watch"]), body.get("name", ""), body.get("rules"), body.get("description", ""), bool(body.get("include_private")),
                                 body.get("manual_include"), body.get("manual_exclude"), body.get("build", ""), body.get("id"))
        except (ValueError, KeyError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.delete("/api/datasets/{dsid}")
    def api_datasets_delete(request: Request, dsid: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        return {"deleted": datasets.delete(dsid)}

    @app.get("/api/datasets/{dsid}/new")
    def api_datasets_new(dsid: int, since: str = ""):
        """Genomes that came into the table after the dataset was last used (or after `since`), and which of them its rules take in."""
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)
        return datasets.newcomers(ds, since)

    @app.post("/api/datasets/{dsid}/add")
    async def api_datasets_add(request: Request, dsid: int):
        """Add genomes to a dataset whatever its rules say (manual include)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        body = await request.json()
        try:
            return datasets.add_manual(dsid, body.get("ids") or [])
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/datasets/{dsid}/qc")
    def api_datasets_qc(dsid: int):
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)
        pub, priv = datasets.member_ids(ds["watch_id"], ds)
        return qc.summary_for(ds["watch_id"], pub, priv)

    @app.post("/api/datasets/{dsid}/qc")
    def api_datasets_qc_run(request: Request, dsid: int, force: int = 0):
        """QC for this dataset's members only (a server job)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)
        return start_job(lambda p: {"redirect": f"/dataset/{dsid}", "result": datasets.run_qc(ds, p, bool(force))}, f"데이터셋 QC · {ds['name']}", "dsqc")

    # ------------------------------------------------------------------------------------------------ preparing private genomes for Pathoplexus (the app never uploads)
    @app.get("/pathoplexus")
    def pathoplexus_page(request: Request, watch: int = 0):
        watches = _rows("SELECT id, name FROM watchlist ORDER BY name")
        wid = watch or next((w["id"] for w in watches if "monkeypox" in w["name"].lower()), watches[0]["id"] if watches else 0)
        return render(request, "pathoplexus.html", nav="pp", watches=watches, wid=wid)

    @app.get("/api/pathoplexus/{watch}/candidates")
    def api_pp_candidates(watch: int, refresh: int = 0):
        if not _rows("SELECT 1 FROM watchlist WHERE id=?", watch):
            return JSONResponse({"error": "unknown virus"}, 404)
        if refresh and pathoplexus.organism_for(watch):
            pathoplexus.template(pathoplexus.organism_for(watch), refresh=True)
        c = pathoplexus.candidates(watch)
        for r in c["rows"]:
            r.pop("file_path", None)
        return c

    @app.post("/api/pathoplexus/{watch}/package")
    async def api_pp_package(request: Request, watch: int):
        """Write the submission files for the chosen private genomes (nothing is uploaded). Body: ids, terms (OPEN|RESTRICTED), restricted_until, group_id, note."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        body = await request.json()
        try:
            return pathoplexus.build_package(watch, [str(i) for i in body.get("ids", [])], body.get("terms", "OPEN"), body.get("restricted_until", ""), str(body.get("group_id") or "").strip(),
                                             str(body.get("note") or ""), bool(body.get("allow_warnings", True)))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/pathoplexus/packages")
    def api_pp_packages(watch: int = 0):
        return pathoplexus.list_packages(watch or None)

    @app.get("/api/pathoplexus/packages/{pid}.zip")
    def api_pp_zip(pid: int):
        try:
            data = pathoplexus.zip_package(pid)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 404)
        return Response(data, media_type="application/zip", headers={"Content-Disposition": f'attachment; filename="pathoplexus_package_{pid}.zip"'})

    @app.post("/api/pathoplexus/packages/{pid}/mark")
    async def api_pp_mark(request: Request, pid: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        body = await request.json()
        try:
            return pathoplexus.mark(pid, body.get("status", ""), body.get("accessions"), body.get("note"))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    # ------------------------------------------------------------------------------------------------ Nextstrain analysis of a dataset
    @app.get("/hub")
    def analysis_hub_page(request: Request, dataset: int = 0, watch: int = 0):
        """The landing page of the analyses: pick a dataset, see how its runs relate to the data today, and open every result view from one place."""
        ds = datasets.get(dataset) if dataset else None
        watches = _rows("SELECT id, name FROM watchlist ORDER BY name")
        wid = ds["watch_id"] if ds else (watch or next((w["id"] for w in watches if "monkeypox" in w["name"].lower()), watches[0]["id"] if watches else 0))
        cur = next((w for w in watches if w["id"] == wid), None)
        return render(request, "hub.html", nav="hub", watches=watches, wid=wid, dsid=ds["id"] if ds else 0, support=support.for_virus(cur["name"]) if cur else [])

    @app.get("/nextstrain")
    def nextstrain_page(request: Request, dataset: int = 0, watch: int = 0):
        ds = datasets.get(dataset) if dataset else None
        watches = _rows("SELECT id, name FROM watchlist ORDER BY name")
        wid = ds["watch_id"] if ds else (watch or next((w["id"] for w in watches if "monkeypox" in w["name"].lower()), watches[0]["id"] if watches else 0))
        return render(request, "nextstrain.html", nav="ns", watches=watches, wid=wid, dsid=ds["id"] if ds else 0)

    @app.get("/api/nextstrain/status")
    def api_ns_status():
        return {**nextstrain.status(), "squirrel": squirrel.status(), "mutprofile": mutprofile.status()}

    @app.post("/api/nextstrain/clone")
    def api_ns_clone(request: Request):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        return start_job(lambda p: {"redirect": "/nextstrain", "result": nextstrain.clone_workflow(p)}, "Nextstrain 워크플로 받기", "nextstrain")

    @app.post("/api/nextstrain/viewer/{action}")
    def api_ns_viewer(request: Request, action: str):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        if action == "stop":
            return {"stopped": nextstrain.stop_viewer()}
        if action != "start":
            return JSONResponse({"error": "start 또는 stop"}, 400)
        return start_job(lambda p: {"redirect": "/nextstrain", "result": {"url": nextstrain.start_viewer(p)}}, "Auspice 보기 서버 켜기", "nextstrain")

    @app.get("/api/datasets/{dsid}/nextstrain")
    def api_ds_nextstrain(dsid: int):
        """Everything the analysis panel needs: the dataset's settings, the effective config, how the kept runs relate to the dataset today."""
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)
        err = ""
        try:
            merged = nextstrain.merged_config(ds)
            import yaml
            merged_text = yaml.safe_dump(merged, allow_unicode=True, sort_keys=False)
            defaults_text = yaml.safe_dump(nextstrain.default_config(ds.get("build")), allow_unicode=True, sort_keys=False)
            psha = nextstrain.params_sha(ds.get("params_yaml", ""))
        except ValueError as e:
            merged_text = defaults_text = psha = ""
            err = str(e)
        now_ids = nextstrain.current_members(ds)
        ccontent = nextstrain.content_sha(ds)
        commit = nextstrain.workflow_commit()
        runs = []
        try:
            ssha, sver = squirrel.settings_sha(ds), squirrel.version() if squirrel.env_prefix() else ""
        except ValueError:
            ssha, sver = "", ""
        try:
            msha_, mver = mutprofile.settings_sha(ds), mutprofile.status()["version"]
        except ValueError:
            msha_, mver = "", ""
        for r in nextstrain.list_runs(dsid):
            if r["tool"] == "mutprofile":
                r["cmp"] = nextstrain.compare(r, now_ids, msha_, mver, ccontent) if msha_ else None
                r["files"] = mutprofile.result_files(r) if r["status"] == "done" else []
                r["summary"] = mutprofile.summary(r) if r["status"] == "done" else {}
                r["url"] = ""
            elif r["tool"] == "clock":
                r["cmp"], r["url"], r["files"] = None, f"/clock/{r['id']}", []
                r["summary"] = {k: v for k, v in clock.summary(r).items() if k in ("n_dated", "clock", "time")} if r["status"] == "done" else {}
            elif r["tool"] == "variantflow":                          # snippy + lofreq of a mutation-profile run: has its own page, is not a Nextstrain run
                r["cmp"], r["url"], r["files"] = None, "", []
            elif r["tool"] == "squirrel":
                try:
                    r["align_only"] = json.loads(r.get("params_yaml") or "{}").get("run_phylo") is False       # the alignment-only run made for the mutation profile: not the Squirrel analysis
                except ValueError:
                    r["align_only"] = False
                r["cmp"] = nextstrain.compare(r, now_ids, ssha, sver, ccontent) if ssha else None
                r["url"] = f"/api/analyses/{r['id']}/files/{r['auspice_file']}" if r["status"] == "done" and r["auspice_file"] else ""
                r["files"] = squirrel.result_files(r) if r["status"] == "done" else []
            else:
                r["cmp"] = nextstrain.compare(r, now_ids, psha, commit, ccontent) if psha else None
                r["url"] = nextstrain.viewer_url(r["name"]) if r["status"] == "done" else ""
            runs.append(r)
        return {"dataset": {"id": ds["id"], "name": ds["name"], "build": ds["build"] or "mpxv", "include_private": ds["include_private"], "params_yaml": ds["params_yaml"], "summary": datasets.describe(ds["rules"])},
                "squirrel": {"settings": squirrel.effective(ds) if ssha else {}, "defaults": squirrel.defaults_for(ds), "saved": ds.get("squirrel_json") or "", "clades": list(squirrel.CLADES)},
                "mutprofile": {"settings": mutprofile.effective(ds) if msha_ else {}, "defaults": mutprofile.defaults_for(ds)},
                "members": len(now_ids), "qc_missing": nextstrain.qc_missing(ds), "members_sha": datasets.members_hash(now_ids), "workflow_commit": commit, "defaults": defaults_text, "merged": merged_text, "error": err, "runs": runs}

    @app.post("/api/datasets/{dsid}/nextstrain/params")
    async def api_ds_ns_params(request: Request, dsid: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        if not datasets.get(dsid):
            return JSONResponse({"error": "unknown dataset"}, 404)
        body = await request.json()
        text = str(body.get("params_yaml") or "")
        try:
            nextstrain.parse_params(text)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        datasets.set_params(dsid, text)
        return {"saved": True}

    @app.post("/api/datasets/{dsid}/squirrel/settings")
    async def api_ds_squirrel_settings(request: Request, dsid: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        if not datasets.get(dsid):
            return JSONResponse({"error": "unknown dataset"}, 404)
        body = await request.json()
        try:
            clean = squirrel.parse_settings(body.get("settings") or {})
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        datasets.set_squirrel(dsid, json.dumps(clean, ensure_ascii=False))
        return {"saved": True, "settings": clean}

    @app.post("/api/datasets/{dsid}/squirrel")
    def api_ds_squirrel_run(request: Request, dsid: int, force: int = 0, cores: int = 0, note: str = ""):
        """Run Squirrel on this dataset (public + private members; a server job; every run is kept)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)

        def work(p):
            r = squirrel.run(ds, cores or None, bool(force), note, p)
            return {"redirect": f"/nextstrain?dataset={dsid}", "result": {"run": r["id"], "name": r["name"], "reused": bool(r.get("reused"))}}
        return start_job(work, f"Squirrel 분석 · {ds['name']}", "squirrel")

    @app.post("/api/datasets/{dsid}/mutprofile/settings")
    async def api_ds_mutprofile_settings(request: Request, dsid: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        if not datasets.get(dsid):
            return JSONResponse({"error": "unknown dataset"}, 404)
        body = await request.json()
        try:
            clean = mutprofile.parse_settings(body.get("settings") or {})
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        datasets.set_mutprofile(dsid, json.dumps(clean, ensure_ascii=False))
        return {"saved": True, "settings": clean}

    @app.post("/api/datasets/{dsid}/mutprofile")
    def api_ds_mutprofile_run(request: Request, dsid: int, force: int = 0, cores: int = 0, note: str = ""):
        """Mutation profile (mutation_profile) with gene / codon annotation for this dataset (public + private; a server job; every run is kept)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)

        def work(p):
            r = mutprofile.run(ds, cores or None, bool(force), note, p)
            return {"redirect": f"/nextstrain?dataset={dsid}", "result": {"run": r["id"], "name": r["name"], "reused": bool(r.get("reused"))}}
        return start_job(work, f"변이 프로파일 · {ds['name']}", "mutprofile")

    @app.get("/vocab")
    def vocab_page(request: Request, watch: int = 0):
        watches = _rows("SELECT id, name FROM watchlist ORDER BY name")
        wid = watch or int(request.cookies.get("gc_watch") or 0) or next((w["id"] for w in watches if "monkeypox" in w["name"].lower()), watches[0]["id"] if watches else 0)
        return render(request, "vocab.html", nav="voc", watches=watches, wid=wid)

    @app.get("/api/vocab/{watch}/preview")
    def api_vocab_preview(watch: int):
        """What the dictionary would do to the curated table as it is now (nothing is written)."""
        w = _rows("SELECT id, name, taxon FROM watchlist WHERE id=?", watch)
        if not w:
            return JSONResponse({"error": "unknown virus"}, 404)
        from .db import connect as _connect
        con = _connect()
        try:
            return {"scope": vocab.scope_for(w[0]["name"], w[0]["taxon"] or ""), "fields": vocab.preview(con, watch, vocab.scope_for(w[0]["name"], w[0]["taxon"] or "")),
                    "dirs": [str(d) for d in vocab.vocab_dirs()], "pending": [dict(r) for r in con.execute("SELECT field, status, COUNT(*) n FROM standardization WHERE watch_id=? GROUP BY field, status", (watch,))]}
        finally:
            con.close()

    def _vocab_scope(watch: int) -> str:
        w = _rows("SELECT name, taxon FROM watchlist WHERE id=?", watch)
        return vocab.scope_for(w[0]["name"], w[0]["taxon"] or "") if w else "*"

    def _vocab_audit(what: str):
        from .db import audit as _audit, connect as _connect
        con = _connect()
        try:
            _audit(con, "vocab", what)
            con.commit()
        finally:
            con.close()

    @app.get("/api/vocab/{watch}/review")
    def api_vocab_review(watch: int):
        """What a person decides: per field the unmapped values (with the closest terms) and the flagged ones."""
        if not _rows("SELECT 1 FROM watchlist WHERE id=?", watch):
            return JSONResponse({"error": "unknown virus"}, 404)
        from .db import connect as _connect
        con = _connect()
        try:
            return {"scope": _vocab_scope(watch), "fields": vocab.review(con, watch, _vocab_scope(watch))}
        finally:
            con.close()

    @app.get("/api/vocab/terms")
    def api_vocab_terms(field: str, q: str = "", watch: int = 0):
        if field not in vocab.FILES:
            return JSONResponse({"error": "unknown field"}, 400)
        return vocab.search_terms(field, q, _vocab_scope(watch) if watch else "*")

    @app.post("/api/vocab/term")
    async def api_vocab_term(request: Request):
        """Add a standard term (optionally with synonyms and a flag) to the user's dictionary."""
        b = await request.json()
        try:
            r = vocab.add_term(b.get("field", ""), b.get("term", ""), b.get("synonyms") or [], b.get("flag", ""), b.get("note", ""), _vocab_scope(int(b.get("watch") or 0)) if b.get("field") == "clade" else "*")
        except (ValueError, KeyError) as e:
            return JSONResponse({"error": str(e)}, 400)
        _vocab_audit(f"용어 추가 {b.get('field')}: {r['term']} (동의어 {len(r['synonyms'])}개)")
        return r

    @app.post("/api/vocab/synonym")
    async def api_vocab_synonym(request: Request):
        b = await request.json()
        try:
            r = vocab.add_synonym(b.get("field", ""), b.get("raw", ""), b.get("term", ""), _vocab_scope(int(b.get("watch") or 0)) if b.get("field") == "clade" else "*")
        except (ValueError, KeyError) as e:
            return JSONResponse({"error": str(e)}, 400)
        _vocab_audit(f"동의어 추가 {b.get('field')}: {b.get('raw')} → {r['term']}")
        return r

    @app.post("/api/vocab/date")
    async def api_vocab_date(request: Request):
        b = await request.json()
        try:
            r = vocab.set_date(b.get("raw", ""), b.get("standard", ""), b.get("note", ""))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        _vocab_audit(f"날짜 지정: {r['raw']} → {r['standard']}")
        return r

    @app.post("/api/vocab/remove")
    async def api_vocab_remove(request: Request):
        b = await request.json()
        if b.get("field") not in vocab.FILES:
            return JSONResponse({"error": "unknown field"}, 400)
        ok = vocab.remove_entry(b["field"], b.get("term", ""), _vocab_scope(int(b.get("watch") or 0)) if b["field"] == "clade" else "*")
        if ok:
            _vocab_audit(f"사용자 사전에서 삭제 {b['field']}: {b.get('term')}")
        return {"removed": ok}

    @app.get("/api/world")
    def api_world():
        """Country outlines (Natural Earth 110m, public domain) for the spread map."""
        return Response(geo.outlines(), media_type="application/json", headers={"Cache-Control": "public, max-age=86400"})

    @app.get("/api/lineages/status")
    def api_lineages_status():
        return lineages.status()

    @app.post("/api/lineages/clone")
    def api_lineages_clone(request: Request):
        """git clone / pull of mpxv-lineages/lineage-designation into data/lineage-designation (a server job)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        return start_job(lambda p: {"result": lineages.clone(p)}, "lineage 정의 받기", "lineages")

    @app.get("/analysis/{rid}/spread")
    def analysis_spread_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "spread.html", nav="res", run=r)

    @app.get("/api/analyses/{rid}/spread")
    def api_analysis_spread(rid: int, view: str = "lineage", bin: str = "auto", outbreak: str = "", min_n: int = 2, lineage: str = "", country: str = "", sites: str = ""):
        """Gantt-like data: lineages (or defining SNPs) x countries x ISO week / month for a mutation-profile run's dataset. `sites` (`7226G,18850C`, positions of the run's reference):
        the carriers of those mutations instead."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return lineages.spread(r, view if view in ("lineage", "mutation") else "lineage", bin, outbreak, max(1, min_n), [x for x in lineage.split(",") if x], [x for x in country.split(",") if x],
                                   lineages.parse_sites(sites) if sites else ())
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/analysis/{rid}/lineage")
    def analysis_lineage_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "lineagetree.html", nav="res", run=r)

    @app.get("/api/analyses/{rid}/lineage_tree")
    def api_analysis_lineage_tree(rid: int, weeks: int = 0):
        """The clade / outbreak / lineage hierarchy with the mutations that introduced each node and the genomes of the run (own, with sub-lineages, within the last `weeks` weeks)."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return lineagetree.tree(r, min(max(0, weeks), 520))
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/analyses/{rid}/lineage_growth")
    def api_analysis_lineage_growth(rid: int, bin: str = "week", weeks: int = 52, country: str = "", min_n: int = 20, ref: str = ""):
        """Growth advantage of the lineages (multinomial logistic regression of their shares over time)."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        from collections import Counter
        try:
            out = lineagetree.growth(r, bin, min(max(4, weeks), 520), country, max(3, min_n), 3, ref)
            out["countries"] = [c for c, _n in Counter(g["country"] for g in lineages._genomes(r) if g["country"]).most_common(40)]
            return out
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/analyses/{rid}/hyphy.zip")
    def api_analysis_hyphy(rid: int, genes: str = ""):
        """Input for HyPhy (FEL / MEME / FUBAR): a codon alignment per requested gene and the tree, public genomes only."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            data = hyphyexport.build(r, [g for g in genes.split(",") if g][:80])
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)
        return Response(data, media_type="application/zip", headers={"Content-Disposition": f'attachment; filename="hyphy-input-run{rid}.zip"'})

    @app.get("/api/analyses/{rid}/clades_tsv")
    def api_analysis_clades_tsv(rid: int, download: int = 0, validate: int = 0):
        """A Nextstrain clades.tsv for the lineages of clade Ib (mpxv-lineages sh2023) in the run's reference coordinates; `validate=1` replays augur's rule on the run's SNV calls and compares with the lineages of the genome table."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            b = cladedefs.build(r)
        except (ValueError, RuntimeError, OSError) as e:
            return JSONResponse({"error": str(e)}, 400)
        if download:
            return Response(b["text"], media_type="text/tab-separated-values", headers={"Content-Disposition": 'attachment; filename="clades.tsv"'})
        out = {k: v for k, v in b.items() if k != "text"}
        out["preview"] = b["text"].split("\n# Lineages")[-1][:3000]
        if validate:
            out["validation"] = cladedefs.validate(r)
        return out

    @app.get("/analysis/{rid}/focus")
    def analysis_focus_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "focus.html", nav="res", run=r)

    @app.get("/api/analyses/{rid}/focus")
    def api_analysis_focus(rid: int, country: str = "South Korea", private: int = 0, max_n: float = -1):
        """One country's genomes (or only the private ones) against the rest of the run: distinct mutations, shared profile, nearest relatives."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return focus.analyze(r, country, bool(private), None if max_n < 0 else min(1.0, max_n))
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/analysis/{rid}/continent")
    def analysis_continent_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "continent.html", nav="res", run=r)

    @app.get("/api/analyses/{rid}/continent")
    def api_analysis_continent(rid: int, by: str = "continent", min_group: int = 10, min_carriers: int = 5, max_n: float = 0.05, adjust: str = "lineage", perms: int = 400, top: int = 80):
        """Mutations by continent (or country): frequencies per group, an overall test, and the same test inside lineages (label permutation), corrected for the number of mutations."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return continents.compare(r, by, max(2, min_group), max(2, min_carriers), min(1.0, max(0.0, max_n)), "lineage" if adjust == "lineage" else "none", min(2000, max(100, perms)), min(300, max(10, top)))
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/analysis/{rid}/concordance")
    def analysis_concordance_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "concordance.html", nav="res", run=r, vc=variantcall.status(), vf=variantflow.status(), nc=nextcladecall.status(r))

    @app.get("/api/analyses/{rid}/concordance")
    def api_analysis_concordance(rid: int, max_n: float = 0.05, source: str = "direct"):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return {**variantcall.concordance_for_run(r, min(1.0, max(0.0, max_n)), source="flow" if source == "flow" else "direct"), "parity": variantflow.parity(r)}

    @app.post("/api/analyses/{rid}/variantcall")
    def api_analysis_variantcall(request: Request, rid: int):
        """snippy + lofreq on the genomes of a finished mutation-profile run (a server job; minutes for a thousand genomes)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        if not variantcall.status()["ok"]:
            return JSONResponse({"error": "snippy / lofreq 환경이 없습니다 (variantcall/*.yml)"}, 400)

        def work(p):
            res = variantcall.run_for_mutprofile(r, progress=lambda i, n: p(100 * i / max(n, 1), f"{i}/{n} 유전체"))
            return {"redirect": f"/analysis/{rid}/concordance", "result": res}
        return start_job(work, f"변이 호출(snippy+lofreq) · {r['name']}", "variantcall")

    @app.post("/api/analyses/{rid}/nextcladecall")
    def api_analysis_nextcladecall(request: Request, rid: int):
        """Nextclade's mutation calls for the genomes of a finished mutation-profile run (a server job; about half a minute per thousand genomes)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        if not nextcladecall.status(r)["ok"]:
            return JSONResponse({"error": "Nextclade 또는 이 바이러스의 Nextclade 데이터셋이 없습니다 (QC 화면에서 데이터셋을 받으세요)"}, 400)

        def work(p):
            return {"redirect": f"/analysis/{rid}/concordance", "result": nextcladecall.run_for_mutprofile(r, progress=p)}
        return start_job(work, f"변이 호출(Nextclade) · {r['name']}", "nextcladecall")

    @app.get("/api/analyses/{rid}/mutation_calls")
    def api_analysis_mutation_calls(rid: int, pos: int, alt: str, max_n: float = 0.05, source: str = "direct"):
        """One mutation tool by tool: how many genomes each caller (the run's own, snippy, lofreq, Nextclade) calls it in, which combinations occur, which tool alone calls it."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return {**variantcall.site_calls(r, pos, alt, min(1.0, max(0.0, max_n)), "flow" if source == "flow" else "direct"),
                "have": {"snippy": variantcall.available(r), "flow": variantcall.available(r, "flow"), "nextclade": nextcladecall.available(r)}}

    @app.get("/api/variantflow/pending")
    def api_variantflow_pending():
        return {"pending": [{"id": r["id"], "name": r["name"], "dataset": r["dataset_name"]} for r in variantflow.pending()], "ok": variantflow.status()["ok"]}

    @app.post("/api/variantflow/all")
    def api_variantflow_all(request: Request):
        """snippy + lofreq (Snakemake flow) for every finished mutation-profile run that has no result yet, as one background job."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        if not variantflow.status()["ok"]:
            return JSONResponse({"error": variantflow.status()["how"]}, 400)
        todo = variantflow.pending()
        if not todo:
            return JSONResponse({"error": "결과가 없는 변이 프로파일 실행이 없습니다 (모두 끝났습니다)."}, 400)
        return start_job(lambda p: {"redirect": "/nextstrain", "result": variantflow.run_all(progress=p)}, f"변이 호출(Snakemake) · 전체 {len(todo)}개 실행", "variantflow")

    @app.post("/api/analyses/{rid}/variantflow")
    def api_analysis_variantflow(request: Request, rid: int):
        """The same snippy + lofreq calls through the Snakemake flow (resumable); kept beside the direct implementation until the flow is complete."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        if not variantflow.status()["ok"]:
            return JSONResponse({"error": variantflow.status()["how"]}, 400)

        def work(p):
            res = variantflow.run_recorded(r, progress=p)
            return {"redirect": f"/analysis/{rid}/concordance", "result": res}
        return start_job(work, f"변이 호출(Snakemake) · {r['name']}", "variantflow")

    @app.get("/analysis/{rid}/insight")
    def analysis_insight_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "insight.html", nav="res", run=r, runref=mutprofile.summary(r).get("reference", ""))

    @app.get("/api/analyses/{rid}/insight")
    def api_analysis_insight(rid: int, window_days: int = 120, min_carriers: int = 3, skip_apobec: int = 0, max_n: float = 0.05):
        """APOBEC3 transmission check, new (non-defining) mutations and clusters, per-lineage signatures of a mutation-profile run."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return insight.insight(r, window_days=max(1, window_days), min_carriers=max(2, min_carriers), skip_apobec=bool(skip_apobec), max_n=min(1.0, max(0.0, max_n)))
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/analysis/{rid}/heatmap")
    def analysis_heatmap_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "heatmap.html", nav="res", run=r)

    @app.get("/api/analyses/{rid}/heatmap")
    def api_analysis_heatmap(rid: int, group_by: str = "lineage", min_n: int = 3, max_n: float = 0.05, min_freq: float = 0.05, min_carriers: int = 2):
        """Lineage (or clade / subcluster) x mutation frequencies with genes, histogram and defining mutations."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return heatmap.build(r, group_by if group_by in ("lineage", "clade", "cluster") else "lineage", max(1, min_n), min(1.0, max(0.0, max_n)), min(1.0, max(0.0, min_freq)), max(1, min_carriers))
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/analyses/{rid}/gvf")
    def api_analysis_gvf(rid: int, group: str = "", group_by: str = "lineage", max_n: float = 0.05):
        """GVF of one lineage (text), or a zip with one GVF per lineage (for VIRUS-MVP / nf-ncov-voc)."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            data = heatmap.gvf(r, group, group_by if group_by in ("lineage", "clade") else "lineage", min(1.0, max(0.0, max_n)))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 404)
        if group:
            return Response(data, media_type="text/plain; charset=utf-8", headers={"Content-Disposition": f'attachment; filename="{group}_annotated.gvf"'})
        return Response(data, media_type="application/zip", headers={"Content-Disposition": 'attachment; filename="gvf_by_lineage.zip"'})

    @app.get("/analysis/{rid}/pca")
    def analysis_pca_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "pca.html", nav="res", run=r)

    @app.get("/api/analyses/{rid}/pca")
    def api_analysis_pca(rid: int, max_n: float = 0.05, min_carriers: int = 3, merge: float = 0.15):
        """Groups of co-occurring mutations (Jaccard clustering) and the principal components of the genomes."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return cooccur.analyse(r, min(1.0, max(0.0, max_n)), max(2, min_carriers), 0.9, min(0.6, max(0.0, merge)))
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/analysis/{rid}/novel")
    def analysis_novel_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "novel.html", nav="res", run=r, runref=mutprofile.summary(r).get("reference", ""))

    @app.get("/api/analyses/{rid}/novel")
    def api_analysis_novel(rid: int, window_days: int = 120, min_size: int = 3, max_n: float = 1.0, homoplasy: int = 0):
        """New-variant discovery: candidate groups sharing mutations no designated lineage is defined by, with status, and the discovery curve of new clade mutations."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return novel.discover(r, max(7, window_days), max(2, min_size), min(1.0, max(0.0, max_n)), bool(homoplasy))
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/analysis/{rid}/outbreak")
    def analysis_outbreak_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "outbreak.html", nav="res", run=r)

    @app.get("/api/analyses/{rid}/outbreak")
    def api_analysis_outbreak(rid: int, min_size: int = 5, hotspot: int = 100, fixed: float = 0.9, max_n: float = 0.05):
        """Genomic subclusters from the tree, cluster-size distribution, APOBEC3 accumulation, per-genome mutation types, gene burden / hotspots / exploratory dN/dS."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return outbreak.analysis(r, max(2, min_size), max(1, hotspot), min(1.0, max(0.5, fixed)), min(1.0, max(0.0, max_n)))

    @app.get("/analysis/{rid}/snipit")
    def analysis_snipit_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return render(request, "snipit.html", nav="res", run=r)

    @app.get("/api/analyses/{rid}/snipit")
    def api_analysis_snipit(rid: int, lineage: str = "", clade: str = "", gene: str = "", start: int = 0, end: int = 0, apobec: int = 0, hide_common: float = 0.0,
                            per_group: int = 8, max_samples: int = 60, group_by: str = "lineage", sort: str = "group"):
        """Data for the snipit-style picture of a mutation-profile run (reference = the run's reference)."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return mutprofile.snipit_data(r, [x for x in lineage.split(",") if x], [x for x in clade.split(",") if x], gene, start, end, bool(apobec), hide_common,
                                          per_group, min(max_samples, 300), group_by if group_by in ("lineage", "clade") else "lineage", sort)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/analyses/{rid}/mutation_gantt")
    def api_analysis_mutation_gantt(rid: int, bin: str = "month", lineage: str = "", clade: str = "", gene: str = "", start: int = 0, end: int = 0, apobec: int = 0, hide_common: float = 0.0,
                                    min_carriers: int = 3, max_rows: int = 150, sort: str = "first", max_n: float = 0.05, defining: str = "all"):
        """The mutations of a run as a Gantt chart: one row per mutation, time in ISO weeks or months, ordered by when each was first observed (collection dates of the carrying genomes)."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        try:
            return mutgantt.build(r, bin, [x for x in lineage.split(",") if x], [x for x in clade.split(",") if x], gene, start, end, bool(apobec), hide_common, max(1, min_carriers),
                                  min(max(1, max_rows), 600), sort if sort in ("first", "pos", "n") else "first", min(1.0, max(0.0, max_n)), defining if defining in ("all", "only", "none") else "all")
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/analyses/{rid}/mutation_series")
    def api_analysis_mutation_series(rid: int, pos: int, alt: str, bin: str = "month", max_n: float = 0.05):
        """One mutation over time: carriers per month/ISO week (by lineage) against the genomes sampled in that period."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "mutprofile" or r["status"] != "done":
            return JSONResponse({"error": "finished mutation-profile run not found"}, 404)
        return mutgantt.series(r, pos, alt, bin, min(1.0, max(0.0, max_n)))

    @app.get("/api/analyses/{rid}/files/{path:path}")
    def api_analysis_file(rid: int, path: str):
        """A file a run left (e.g. Squirrel's HTML report, tree, alignment), only from inside that run's folder."""
        r = nextstrain.get_run(rid)
        if not r or not r["run_dir"]:
            return JSONResponse({"error": "unknown run"}, 404)
        base = Path(r["run_dir"]).resolve()
        f = (base / path).resolve()
        if base not in f.parents or not f.is_file():
            return JSONResponse({"error": "no such file"}, 404)
        return FileResponse(f, filename=None if f.suffix in (".html", ".csv", ".tsv", ".txt", ".json") else f.name)

    @app.post("/api/datasets/{dsid}/analyze")
    def api_ds_analyze(request: Request, dsid: int, force: int = 0, runtime: str = "ambient", cores: int = 0, note: str = ""):
        """Run the Nextstrain mpox workflow on this dataset (a server job; every run is kept)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)

        def work(p):
            r = nextstrain.run(ds, runtime, cores or None, bool(force), note, p)
            return {"redirect": f"/nextstrain?dataset={dsid}", "result": {"run": r["id"], "name": r["name"], "reused": bool(r.get("reused"))}}
        return start_job(work, f"Nextstrain 분석 · {ds['name']}", "nextstrain")

    @app.post("/api/datasets/{dsid}/analyze_all")
    def api_ds_analyze_all(request: Request, dsid: int, steps: str = "", force: int = 0, cores: int = 0):
        """QC, Squirrel, mutation profile, variant calls and Nextstrain for this dataset, one after the other, as one background job (steps unchanged since the last run are reused)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)
        want = [s for s in steps.split(",") if s in pipeline.STEPS] or None
        todo, _notes = pipeline.plan(ds, want)
        if not todo:
            return JSONResponse({"error": "실행할 단계가 없습니다."}, 400)
        return start_job(lambda p: {"redirect": f"/dataset/{dsid}", "result": pipeline.run_all(ds, want, bool(force), cores or None, p)}, f"전체 분석 · {ds['name']} ({len(todo)}단계)", "pipeline")

    @app.post("/api/datasets/{dsid}/clock")
    async def api_ds_clock(request: Request, dsid: int, force: int = 0):
        """Molecular clock and tMRCA (root-to-tip regression + time-scaled tree) on this dataset's tree (a server job). Settings in the JSON body (all optional)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)
        try:
            body = await request.json()
        except ValueError:
            body = {}
        try:
            s = clock.parse_settings(body.get("settings") if isinstance(body, dict) else None)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

        def work(p):
            r = clock.run(ds, s, bool(force), "", p)
            return {"redirect": f"/clock/{r['id']}", "result": {"run": r["id"], "name": r["name"], "reused": bool(r.get("reused"))}}
        return start_job(work, f"분자시계·tMRCA · {ds['name']}", "clock")

    @app.get("/clock/{rid}")
    def clock_page(request: Request, rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        return render(request, "clock.html", nav="res", run=r)

    @app.get("/api/clock/{rid}/beast")
    def api_clock_beast(rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        b = [x for x in nextstrain.list_runs(r["dataset_id"]) if x["tool"] == "beast" and str(x.get("note") or "") == f"clock:{rid}"]
        d = beastprep.folder(r)
        return {"status": beastprep.status(), "defaults": beastprep.DEFAULTS, "runs": [{k: x[k] for k in ("id", "name", "status", "started", "finished", "error")} for x in b],
                "result": beastprep.result(r), "files": [{"name": f.name, "size": f.stat().st_size} for f in sorted(d.iterdir()) if f.is_file()] if d.exists() else [],
                "tips_usable": sum(1 for t in clock.summary(r).get("tips", []) if not t["outlier"] and t["precision"] == "day"),
                "trees": sorted(f.parent.name for f in d.glob("*/beast.trees")) + (["" ] if (d / "beast.trees").exists() else [])}

    @app.get("/api/clock/{rid}/beast/tree")
    def api_clock_beast_tree(rid: int, tag: str = "", kind: str = "hipstr", file: str = "tree", download: int = 0):
        """The summary tree of a BEAST chain (NEXUS), or `file=tips` for the tip table; made on first use with TreeAnnotator."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        try:
            info = beastprep.mcc_tree(r, tag, kind)
            if file == "annotated":
                return FileResponse(beastprep.annotated_tree(r, tag, kind), media_type="text/plain; charset=utf-8")
        except (ValueError, RuntimeError) as e:
            return JSONResponse({"error": str(e)}, 400)
        name = info["tips"] if file == "tips" else info["tree"]
        if not name:
            return JSONResponse({"error": "not available"}, 404)
        path = beastprep.folder(r, tag) / name
        headers = {"X-Newest-Tip": str(info["newest_tip"] or "")}
        return FileResponse(path, media_type="text/plain; charset=utf-8", filename=(f"{r['name']}-{tag or 'beast'}-{name}" if download else None), headers=headers)

    @app.get("/beasttree/{rid}")
    def beast_tree_page(request: Request, rid: int, tag: str = ""):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        return render(request, "beasttree.html", nav="res", run=r, tag=tag, pt=peartree.status())

    @app.post("/api/peartree/install")
    def api_peartree_install():
        return start_job(lambda p: peartree.install(p), "PearTree 설치", "install")

    @app.get("/tools/peartree/{path:path}")
    def peartree_files(path: str):
        base = peartree.root().resolve()
        f = (base / path).resolve()
        if not peartree.installed() or base not in f.parents or not f.is_file():
            return JSONResponse({"error": "not found"}, 404)
        return FileResponse(f)

    @app.post("/api/clock/{rid}/beast/prepare")
    async def api_clock_beast_prepare(request: Request, rid: int):
        """Write the BEAST package (subsample, alignment, dates, XML) without running BEAST."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        try:
            body = await request.json()
        except ValueError:
            body = {}
        try:
            return beastprep.prepare(r, body.get("settings") if isinstance(body, dict) else None)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.post("/api/clock/{rid}/beast/run")
    async def api_clock_beast_run(request: Request, rid: int):
        """Prepare the package and run BEAST on it (a long background job: hours for a long chain)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        if not beastprep.status()["ok"]:
            return JSONResponse({"error": beastprep.status()["how"]}, 400)
        try:
            body = await request.json()
        except ValueError:
            body = {}
        try:
            s = beastprep.parse_settings(body.get("settings") if isinstance(body, dict) else None)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        return start_job(lambda p: {"redirect": f"/clock/{rid}", "result": beastprep.run_recorded(r, s, p)}, f"BEAST · {r['name']}", "beast")

    @app.get("/api/clock/{rid}/bets")
    def api_clock_bets(rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        b = [x for x in nextstrain.list_runs(r["dataset_id"]) if x["tool"] == "bets" and str(x.get("note") or "") == f"clock:{rid}"]
        return {"status": beastprep.status(), "defaults": beastprep.BETS_DEFAULTS, "runs": [{k: x[k] for k in ("id", "name", "status", "started", "finished", "error")} for x in b], "result": beastprep.bets_result(r)}

    @app.post("/api/clock/{rid}/bets/run")
    async def api_clock_bets_run(request: Request, rid: int):
        """BETS (is there a temporal signal?): two BEAST analyses of the same model, with the real dates and with all genomes at one date, and the Bayes factor of their marginal likelihoods (a long background job)."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        if not beastprep.status()["ok"]:
            return JSONResponse({"error": beastprep.status()["how"]}, 400)
        try:
            body = await request.json()
        except ValueError:
            body = {}
        try:
            s = beastprep.parse_bets_settings(body.get("settings") if isinstance(body, dict) else None)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        return start_job(lambda p: {"redirect": f"/clock/{rid}", "result": beastprep.run_bets_recorded(r, s, p)}, f"BETS · {r['name']}", "bets")

    @app.get("/api/clock/{rid}/rtt_class")
    def api_clock_rtt_class(rid: int):
        """Root-to-tip regression of the dataset's Nextstrain tree once for APOBEC3-type and once for the other substitutions (and per lineage)."""
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        try:
            return clockclass.analyze(r)
        except (ValueError, RuntimeError, KeyError) as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/clock/{rid}/models")
    def api_clock_models(rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        b = [x for x in nextstrain.list_runs(r["dataset_id"]) if x["tool"] == "bets" and str(x.get("note") or "") == f"clock:{rid}:models"]
        return {"status": beastprep.status(), "defaults": beastprep.BETS_DEFAULTS, "runs": [{k: x[k] for k in ("id", "name", "status", "started", "finished", "error")} for x in b], "result": beastprep.bets_result(r, True)}

    @app.post("/api/clock/{rid}/models/run")
    async def api_clock_models_run(request: Request, rid: int):
        """Compare the four clock x tree-prior models (strict / relaxed x constant / exponential) by their marginal likelihoods (GSS and / or PS/SS), all with the real dates: a long background job."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        if not beastprep.status()["ok"]:
            return JSONResponse({"error": beastprep.status()["how"]}, 400)
        try:
            body = await request.json()
        except ValueError:
            body = {}
        try:
            s = beastprep.parse_bets_settings({**((body.get("settings") if isinstance(body, dict) else None) or {}), "compare": True})
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        return start_job(lambda p: {"redirect": f"/clock/{rid}", "result": beastprep.run_bets_recorded(r, s, p)}, f"모델 비교 · {r['name']}", "bets")

    @app.get("/api/clock/{rid}")
    def api_clock(rid: int):
        r = nextstrain.get_run(rid)
        if not r or r["tool"] != "clock":
            return JSONResponse({"error": "clock run not found"}, 404)
        return {"run": {k: r[k] for k in ("id", "name", "status", "started", "finished", "error", "dataset_id")}, "summary": clock.summary(r), "files": clock.result_files(r)}

    @app.get("/api/analyses")
    def api_analyses(dataset: int = 0, watch: int = 0):
        out = nextstrain.list_runs(dataset or None, watch or None)
        for r in out:
            r["url"] = nextstrain.viewer_url(r["name"]) if r["status"] == "done" else ""
        return out

    @app.get("/api/analyses/{rid}")
    def api_analysis(rid: int, tail: int = 60):
        r = nextstrain.get_run(rid)
        if not r:
            return JSONResponse({"error": "unknown run"}, 404)
        lp = Path(r["log_path"] or "")
        r["log_tail"] = "".join(lp.read_text(encoding="utf-8", errors="replace").splitlines(True)[-max(1, min(tail, 500)):]) if lp.exists() else ""
        r["url"] = nextstrain.viewer_url(r["name"]) if r["status"] == "done" else ""
        return r

    @app.post("/api/analyses/{rid}/note")
    async def api_analysis_note(request: Request, rid: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        body = await request.json()
        nextstrain.set_note(rid, str(body.get("note") or ""))
        return {"saved": True}

    @app.delete("/api/analyses/{rid}")
    def api_analysis_delete(request: Request, rid: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        return {"deleted": nextstrain.delete_run(rid)}

    @app.get("/api/datasets/{dsid}/history")
    def api_datasets_history(dsid: int):
        return datasets.export_history(dsid)

    @app.post("/api/datasets/{dsid}/export")
    def api_datasets_export(request: Request, dsid: int):
        """Write this dataset's files (metadata.tsv, sequences.fasta, manifest, config; plus the private files if the dataset includes them) as a server job."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        ds = datasets.get(dsid)
        if not ds:
            return JSONResponse({"error": "unknown dataset"}, 404)
        return start_job(lambda p: {"redirect": f"/datasets?watch={ds['watch_id']}", "result": export.write(ds["watch_id"], None, "nextstrain", p, dataset=ds)},
                         f"데이터셋 내보내기 · {ds['name']}", "dataset")

    @app.get("/api/export/profiles")
    def api_export_profiles():
        return {k: {"label": v["label"], "id_column": v["id"], "columns": v["columns"]} for k, v in export.PROFILES.items()}

    @app.get("/api/export/{watch}/preview")
    def api_export_preview(request: Request, watch: int, profile: str = "nextstrain"):
        if not _export_watch(watch):
            return JSONResponse({"error": "unknown virus"}, 404)
        try:
            return export.preview(watch, profile, include_private=request.query_params.get("include_private") in ("1", "true"), **_xf(request))
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/export/{watch}/metadata.tsv")
    def api_export_metadata(request: Request, watch: int, profile: str = "nextstrain"):
        """The metadata table as TSV - the same rows, in the same order, as sequences.fasta."""
        w = _export_watch(watch)
        if not w:
            return JSONResponse({"error": "unknown virus"}, 404)
        try:
            ds = _ds(request, watch)
            cols, rows, _, _ = export.table(watch, profile, ds, **({} if ds else _xf(request)))
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
            ds = _ds(request, watch)
            _, _, sel, _ = export.table(watch, profile, ds, **({} if ds else _xf(request)))
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
        inc = request.query_params.get("include_private") in ("1", "true")
        return start_job(lambda p: {"redirect": "/export", "result": export.write(watch, None, profile, p, include_private=inc, **flt)},
                         f"내보내기 · {w['name']} ({profile})" + (" + 비공개" if inc else ""), "export")

    # ------------------------------------------------------------------------------------------------ private store
    @app.get("/private")
    def private_page(request: Request, watch: int = 0):
        watches = _rows("SELECT id, name FROM watchlist ORDER BY name")
        wid = watch or (watches[0]["id"] if watches else 0)
        tx = _rows("SELECT taxon FROM watchlist WHERE id=?", wid)
        return render(request, "private.html", nav="priv", watches=watches, wid=wid, rows=private.list_rows(wid) if wid else [],
                      hist=private.history(wid) if wid else [], taxon=(tx[0]["taxon"] if tx else ""), clades=private.CLADE_CHOICES.get((tx[0]["taxon"] if tx else "").lower(), []),
                      columns=private.COLUMNS, min_cov=config.MIN_COVERAGE_PCT)

    async def _private_form(request: Request):
        form = await request.form()
        return form

    async def _private_files(request: Request):
        """-> (metadata bytes, fasta text, filename, defaults, note) from the multipart form; metadata may be .xlsx / .tsv / .csv."""
        form = await request.form()
        m, f = form.get("metadata"), form.get("fasta")
        if m is None or f is None or not hasattr(m, "read"):
            return None, None, "", {}, ""
        defaults = {"clade": str(form.get("default_clade") or "").strip()}
        return await m.read(), (await f.read()).decode("utf-8", "replace"), getattr(m, "filename", "") or "", defaults, str(form.get("note") or "").strip()

    @app.post("/api/private/parse")
    async def api_private_parse(request: Request):
        """Read an uploaded spreadsheet (.xlsx) / TSV / CSV into columns and rows, so it can be opened in the page's sheet and fixed there. Nothing is stored."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        form = await request.form()
        m = form.get("metadata")
        if m is None or not hasattr(m, "read"):
            return JSONResponse({"error": "파일이 필요합니다."}, 400)
        try:
            cols, rows = private._read_table(await m.read(), getattr(m, "filename", "") or "")
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)
        return {"cols": cols, "rows": rows[:5000]}

    @app.get("/api/private/template.xlsx")
    def api_private_template(watch: int = 0):
        w = _rows("SELECT name, taxon FROM watchlist WHERE id=?", watch)
        try:
            data = private.template_xlsx(w[0]["name"] if w else None, w[0]["taxon"] if w else None)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 500)
        return Response(data, media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        headers={"Content-Disposition": 'attachment; filename="private_genomes_template.xlsx"'})

    @app.get("/api/private/template.tsv")
    def api_private_template_tsv():
        return Response(private.template_tsv(), media_type="text/tab-separated-values; charset=utf-8",
                        headers={"Content-Disposition": 'attachment; filename="private_genomes_template.tsv"'})

    @app.post("/api/private/{watch}/check")
    async def api_private_check(request: Request, watch: int):
        """Dry run of an upload: multipart `metadata` (.xlsx / .tsv / .csv) and `fasta`, optional `default_clade`. Nothing is stored."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        mt, fa, name, dflt, _ = await _private_files(request)
        if mt is None:
            return JSONResponse({"error": "metadata와 fasta 파일이 필요합니다."}, 400)
        try:
            rep = private.check(watch, mt, fa, None, dflt, name)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 404 if "찾을 수 없" in str(e) else 400)
        return {k: rep[k] for k in ("ok", "n", "new", "update", "unchanged", "id_column", "errors", "warnings")}

    @app.post("/api/private/{watch}/import")
    async def api_private_import(request: Request, watch: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        mt, fa, name, dflt, note = await _private_files(request)
        if mt is None:
            return JSONResponse({"error": "metadata와 fasta 파일이 필요합니다."}, 400)
        try:
            return private.import_texts(watch, mt, fa, name, None, dflt, note)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.post("/api/private/{watch}/one")
    async def api_private_one(request: Request, watch: int, dry_run: int = 0):
        """One genome: form fields (PPX_accession, strain, date, country, ...) plus `sequence` (text) or `fasta` (file). dry_run=1 only checks."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        form = await request.form()
        seq = str(form.get("sequence") or "")
        f = form.get("fasta")
        if f is not None and hasattr(f, "read"):
            seq = (await f.read()).decode("utf-8", "replace") or seq
        skip = {"sequence", "fasta", "default_clade", "note"}
        fields = {k: str(v) for k, v in form.items() if k not in skip and not hasattr(v, "read")}
        try:
            return private.one(watch, fields, seq, "(하나씩 입력)", bool(dry_run), {"clade": str(form.get("default_clade") or "").strip()}, str(form.get("note") or "").strip())
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.get("/api/private/{watch}")
    def api_private_list(watch: int):
        return private.list_rows(watch)

    @app.get("/api/private/{watch}/history")
    def api_private_history(watch: int, id: str = ""):
        return private.history(watch, id or None)

    @app.post("/api/private/{watch}/{pid}/restore")
    def api_private_restore(request: Request, watch: int, pid: str, version: int):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        try:
            return private.restore(watch, pid, version)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.delete("/api/private/{watch}/{pid}")
    def api_private_delete(request: Request, watch: int, pid: str):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        return {"deleted": private.delete(watch, pid)}

    @app.post("/api/private/{watch}/qc")
    def api_private_qc(request: Request, watch: int, force: int = 0):
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        w = _rows("SELECT name FROM watchlist WHERE id=?", watch)
        if not w:
            return JSONResponse({"error": "unknown virus"}, 404)
        return start_job(lambda p: {"redirect": f"/private?watch={watch}", "result": qc.run_private(watch, p, bool(force))}, f"비공개 QC · {w[0]['name']}", "private")

    @app.get("/api/exports")
    def api_exports(virus: str = ""):
        return export.list_exports(virus or None)

    # ------------------------------------------------------------------------------------------------ QC (sequence stats, Nextclade)
    @app.get("/qc")
    def qc_page(request: Request, watch: int = 0):
        watches = _rows("SELECT w.id, w.name, w.taxon, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id AND IFNULL(g.genome_status,'')!='Removed') n, "
                       "(SELECT COUNT(*) FROM genome_sequence s JOIN genome g ON g.genome_id=s.genome_id WHERE g.watch_id=w.id) ns FROM watchlist w ORDER BY w.name")
        return render(request, "qc.html", nav="qc", watches=watches, tool_version=qc.nextclade_version(), tool_ok=qc.nextclade_bin() is not None, checkv=qc.checkv_status(),
                      summary={w["id"]: qc.summary(w["id"]) for w in watches},
                      datasets_map={w["id"]: qc.dataset_for(w["taxon"]) for w in watches})

    @app.get("/api/qc/{watch}/squirrel")
    def api_qc_squirrel(watch: int):
        """Squirrel's sequence QC of the latest finished Squirrel run, next to Nextclade's."""
        return qc.squirrel_qc(watch)

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
        return render(request, "analysis.html", nav="sup", st=analysis.status())

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

    # ---- cluster: the page and its administration (only from this machine unless GC_CLUSTER_ADMIN_REMOTE=1), and the API the workers use (a key per worker)
    def _cluster_admin(request: Request):
        """-> an error response, or None. Creating keys and queueing work is the owner's business: refused from other machines by default."""
        if not _same_origin(request):
            return JSONResponse({"error": "cross-site request refused"}, 403)
        host = request.client.host if request.client else ""
        if host not in ("127.0.0.1", "::1", "localhost", "testclient") and os.environ.get("GC_CLUSTER_ADMIN_REMOTE") != "1":
            return JSONResponse({"error": "the cluster can be administered from the main server only"}, 403)
        return None

    def _worker(request: Request):
        h = request.headers.get("authorization", "")
        return cluster.authenticate(h[7:] if h.lower().startswith("bearer ") else None)

    @app.get("/cluster")
    def cluster_page(request: Request):
        return render(request, "cluster.html", nav="cluster", lease=cluster.LEASE)

    @app.get("/api/cluster/state")
    def cluster_state():
        return cluster.state()

    @app.get("/api/cluster/job/{jid}")
    def cluster_job_get(jid: str):
        j = cluster.get_job(jid)
        return j or JSONResponse({"error": "no such job"}, 404)

    @app.post("/api/cluster/nodes")
    async def cluster_node_create(request: Request):
        if (e := _cluster_admin(request)):
            return e
        b = await request.json()
        try:
            return cluster.create_node(str(b.get("name", "")), bool(b.get("allow_private")), str(b.get("note", "")))
        except ValueError as ex:
            return JSONResponse({"error": str(ex)}, 400)

    @app.post("/api/cluster/nodes/{nid}/revoke")
    def cluster_node_revoke(request: Request, nid: str):
        if (e := _cluster_admin(request)):
            return e
        return {"revoked": cluster.revoke_node(nid)}

    @app.post("/api/cluster/jobs")
    async def cluster_job_submit(request: Request):
        if (e := _cluster_admin(request)):
            return e
        b = await request.json()
        try:
            return {"id": cluster.submit_beast(config.DATA_DIR / str(b["folder"]) if not Path(str(b["folder"])).is_absolute() else b["folder"], str(b.get("label", "")),
                                               int(b.get("max_attempts", 3)))}
        except (ValueError, KeyError) as ex:
            return JSONResponse({"error": str(ex)}, 400)

    @app.post("/api/cluster/jobs/{jid}/cancel")
    def cluster_job_cancel(request: Request, jid: str):
        if (e := _cluster_admin(request)):
            return e
        return {"cancelled": cluster.cancel(jid)}

    @app.post("/api/cluster/w/hello")
    async def cluster_w_hello(request: Request):
        n = _worker(request)
        if not n:
            return JSONResponse({"error": "unknown worker key"}, 403)
        return cluster.hello(n, await request.json())

    @app.post("/api/cluster/w/claim")
    async def cluster_w_claim(request: Request):
        n = _worker(request)
        if not n:
            return JSONResponse({"error": "unknown worker key"}, 403)
        return {"job": cluster.claim(n, (await request.json()).get("kinds"))}

    @app.post("/api/cluster/w/jobs/{jid}/heartbeat")
    async def cluster_w_beat(request: Request, jid: str):
        n = _worker(request)
        if not n:
            return JSONResponse({"error": "unknown worker key"}, 403)
        b = await request.json()
        return cluster.heartbeat(n, jid, b.get("pct"), b.get("msg"))

    @app.get("/api/cluster/w/jobs/{jid}/in/{name}")
    def cluster_w_input(request: Request, jid: str, name: str):
        n = _worker(request)
        if not n:
            return JSONResponse({"error": "unknown worker key"}, 403)
        try:
            return FileResponse(cluster.input_path(n, jid, name))
        except (PermissionError, FileNotFoundError, ValueError) as ex:
            return JSONResponse({"error": str(ex)}, 409)

    @app.put("/api/cluster/w/jobs/{jid}/out/{name}")
    async def cluster_w_output(request: Request, jid: str, name: str, sha: str = ""):
        n = _worker(request)
        if not n:
            return JSONResponse({"error": "unknown worker key"}, 403)
        tmp = cluster.root() / "tmp" / f"{jid}-{os.getpid()}-{time.time_ns()}"
        tmp.parent.mkdir(parents=True, exist_ok=True)
        with open(tmp, "wb") as f:
            async for chunk in request.stream():
                f.write(chunk)
        try:
            return cluster.store_output(n, jid, name, tmp, sha)
        except PermissionError as ex:
            return JSONResponse({"error": str(ex)}, 409)
        except ValueError as ex:
            return JSONResponse({"error": str(ex)}, 400)

    @app.post("/api/cluster/w/jobs/{jid}/done")
    async def cluster_w_done(request: Request, jid: str):
        n = _worker(request)
        if not n:
            return JSONResponse({"error": "unknown worker key"}, 403)
        b = await request.json()
        try:
            return cluster.complete(n, jid, b.get("outputs", []), b.get("seconds"))
        except PermissionError as ex:
            return JSONResponse({"error": str(ex)}, 409)
        except ValueError as ex:
            return JSONResponse({"error": str(ex)}, 400)

    @app.post("/api/cluster/w/jobs/{jid}/fail")
    async def cluster_w_fail(request: Request, jid: str):
        n = _worker(request)
        if not n:
            return JSONResponse({"error": "unknown worker key"}, 403)
        return cluster.fail(n, jid, str((await request.json()).get("error", "")))

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
            msg = f"서열 {res['fetched']}건 다운로드 완료 (다음: QC 탭에서 QC 실행)"
            if res["failed"]:
                first = next(iter(res["failed"].items()))
                msg += f", {len(res['failed'])}건 실패 (예: {first[0]} - {first[1]})"
            return f"/genomes?watch={watch}&msg={msg}"
        return start_job(work, "서열(FASTA) 받기", "fasta")

    @app.post("/changeset/{cid}/dryrun")
    def rerun(cid: int):
        curate.dryrun(cid)
        return RedirectResponse(f"/changeset/{cid}", 303)

    def _genome_view(watch: int, accessions: set | None = None, private_ids: set | None = None, cols=None, facets=None):
        """The rows and settings the genome table page needs. accessions / private_ids narrow it to a dataset's members (private rows are
        added, marked as private); both None = every curated genome of the virus."""
        watches = _rows("SELECT w.id, w.name, w.taxon_id, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id AND IFNULL(g.genome_status,'')!='Removed') n FROM watchlist w")
        cur = next((w for w in watches if w["id"] == watch), None)
        rows = _rows("SELECT g.*, s.md5, s.gc_content FROM genome g LEFT JOIN genome_sequence s ON s.accession=g.accession "
                     "WHERE g.watch_id=? AND IFNULL(g.genome_status,'')!='Removed' ORDER BY g.genome_id", watch)
        if accessions is not None:
            rows = [r for r in rows if r["accession"] in accessions]
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
        exm = exclusions.marked(watch, "public")
        for n, r in enumerate(rows):
            q = qcm.get(r["accession"], {})
            r["excluded"] = "제외" if r["accession"] in exm else "포함"
            r["excl_reason"] = (exm.get(r["accession"]) or {}).get("reason") or ""
            r["qc"], r["nc_clade"], r["nc_lineage"] = q.get("qc", "미실행"), q.get("nc_clade"), q.get("nc_lineage")
            r.pop("taxon_lineage_ids", None)
            r.pop("taxon_lineage_names", None)
            r["id"] = n
            r["genbank"] = f"{r['accession']}.{r['accession_version']}" if r.get("accession_version") else r["accession"]
            r["_prov"] = prov.get(r["genome_id"], {})
            r["has_sra"] = "Yes" if r.get("sra_accession") else "No"
            r["has_fasta"] = "Yes" if r.get("md5") else "No"
            r["completeness"], r["ref_coverage"] = completeness.classify(r.get("genome_status"), r.get("genome_length"), ref and ref[1])
            r["source_kind"] = "공개"
        if private_ids:
            for rec in private.rows_for_export(watch):
                if rec["id"] in private_ids:
                    pr = export.private_as_genome(rec)
                    pr["id"] = len(rows)
                    rows.append(pr)
        ppx = _ppx_slug(watch, [r for r in rows if not r.get("is_private")])
        data = {"rows": rows, "name": "genomes", "watch": cur and cur["name"], "lineage": lineage, "ppx": ppx,
                "ref_note": completeness.describe(ref),
                "cols": _cols(cols or GN_COLS, rows), "facets": _facets(facets or GN_FACETS, rows)}
        return watches, data

    @app.get("/api/watch/{wid}/intake")
    def api_watch_intake(wid: int, bin: str = "day", cbin: str = "month", sel: str = ""):
        """Intake history of a virus's genomes (per day / week / month) and their collection dates, optionally only of the genomes that came in during one intake period."""
        return intake.history(wid, bin, cbin, sel)

    @app.get("/genomes")
    def genomes(request: Request, watch: int = 0, msg: str = ""):
        watch = current_watch(request, watch, _rows("SELECT w.id, w.name, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n FROM watchlist w"))
        # the private genomes of this virus are part of the table (marked in the *출처* column and filterable), like in a dataset that includes them
        priv = {rec["id"] for rec in private.rows_for_export(watch)} if watch else set()
        cols = GN_COLS[:2] + [("source_kind", "출처", True, True)] + GN_COLS[2:]
        facets = [("source_kind", "출처 (공개/비공개)")] + GN_FACETS
        watches, data = _genome_view(watch, None, priv, cols, facets)
        return render(request, "genomes.html", data=data, nav="gen", watches=watches, watch=watch, msg=msg, private_n=len(priv), dataset=None)

    @app.get("/dataset/{dsid}")
    def dataset_view(request: Request, dsid: int, msg: str = ""):
        """A saved dataset browsed like the genome table: filters, sorting, detail panel, TSV; private members (if the dataset includes them) are marked."""
        ds = datasets.get(dsid)
        if not ds:
            return RedirectResponse("/datasets", 303)
        pub, priv = datasets.member_ids(ds["watch_id"], ds)
        cols = GN_COLS[:2] + [("source_kind", "출처", True, True)] + GN_COLS[2:]
        facets = [("source_kind", "출처 (공개/비공개)")] + GN_FACETS
        watches, data = _genome_view(ds["watch_id"], pub, priv, cols, facets)
        data["name"] = "dataset_" + export.slug(ds["name"])
        n_seq = sum(1 for r in data["rows"] if r.get("has_fasta") == "Yes")
        info = {"public": len(pub), "private": len(priv), "no_seq": len(pub) - sum(1 for r in data["rows"] if not r.get("is_private") and r.get("has_fasta") == "Yes"), "with_seq": n_seq}
        return render(request, "genomes.html", data=data, nav="ds", watches=watches, watch=ds["watch_id"], msg=msg, private_n=0, dataset=ds, ds_info=info,
                      ds_summary=datasets.describe(ds["rules"]), ds_qc=qc.summary_for(ds["watch_id"], pub, priv))

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
                          "LEFT JOIN genome_sequence s ON s.accession=g.accession WHERE IFNULL(g.genome_status,'')!='Removed' ORDER BY g.genome_id"), "genomes.tsv")

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
    def dashboard(request: Request, watch: int = 0, dataset: int = 0):
        watches = _rows("SELECT w.id, w.name, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n FROM watchlist w")
        ds = datasets.get(dataset) if dataset else None
        watch = ds["watch_id"] if ds else current_watch(request, watch, watches)
        g = _rows("SELECT * FROM genome WHERE watch_id=? AND IFNULL(genome_status,'')!='Removed'", watch)
        pub = priv = None
        if ds:                                   # the same charts, over the dataset's members only (private ones too when the dataset includes them)
            pub, priv = datasets.member_ids(watch, ds)
            g = [r for r in g if r["accession"] in pub]
        _c = connect()
        try:
            _ref = completeness.reference_map(_c).get(watch)
        finally:
            _c.close()
        for r in g:
            r["completeness"], _ = completeness.classify(r.get("genome_status"), r.get("genome_length"), _ref and _ref[1])
        if ds and priv:
            g += [export.private_as_genome(rec) for rec in private.rows_for_export(watch) if rec["id"] in priv]
        prov = _rows("SELECT p.source, COUNT(*) n FROM field_provenance p JOIN genome g ON g.genome_id=p.genome_id "
                     "WHERE g.watch_id=? GROUP BY p.source", watch) if not ds else \
            _rows("SELECT p.source, g.accession FROM field_provenance p JOIN genome g ON g.genome_id=p.genome_id WHERE g.watch_id=?", watch)
        if ds:
            from collections import Counter
            cnt = Counter(r["source"] for r in prov if r["accession"] in pub)
            prov = [{"source": k, "n": v} for k, v in cnt.items()]
        last = _rows("SELECT last_checked FROM watchlist WHERE id=?", watch)
        return render(request, "dashboard.html", nav="dash", watches=watches, watch=watch, dataset=ds,
                      datasets=datasets.listing(watch), data=dashboard_data(g, prov, last[0]["last_checked"] if last else None))

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
