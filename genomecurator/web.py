import json
import threading
import uuid
from pathlib import Path

from fastapi import FastAPI, Form, Request
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse, Response
from fastapi.templating import Jinja2Templates

from . import config, curate
from .db import connect

T = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))
JOBS: dict[str, dict] = {}


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


def start_job(fn) -> JSONResponse:
    """Run fn(progress) -> redirect_url in a thread; the page polls /api/job/<id> and shows a toast."""
    jid = uuid.uuid4().hex[:8]
    JOBS[jid] = {"pct": 0, "msg": "시작 중…", "done": False, "error": None, "redirect": "/"}

    def prog(pct, msg):
        JOBS[jid].update(pct=pct, msg=msg)

    def run():
        try:
            JOBS[jid]["redirect"] = fn(prog)
            JOBS[jid]["pct"] = 100
        except BaseException as e:
            JOBS[jid]["error"] = str(e) or e.__class__.__name__
        JOBS[jid]["done"] = True

    threading.Thread(target=run, daemon=True).start()
    return JSONResponse({"job": jid})


CS_COLS = [("kind", "종류", True, True), ("decision", "결정", True, True), ("integrity", "정합성", True, True),
           ("conflict_fields", "불일치 필드", True, False), ("accession", "Accession", True, False),
           ("accession_version", "Ver", True, False), ("genome_name", "Genome Name", True, False),
           ("genome_status", "Status", True, False), ("genome_length", "Length", True, False),
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
CS_FACETS = [("kind", "종류"), ("decision", "결정"), ("integrity", "정합성(소스 간 일치)"), ("genome_status", "Genome Status"), ("host_name", "Host"),
             ("isolation_country", "Country"), ("collection_year", "Collection Year"), ("has_sra", "SRA")]
GN_COLS = [("genome_name", "Genome Name", True), ("strain", "Strain", True), ("genbank", "GenBank Accessions", True),
           ("genome_length", "Size", True), ("genome_status", "Genome Status", True), ("collection_date", "Collection Date", True),
           ("isolation_country", "Isolation Country", True), ("host_name", "Host Name", True),
           ("clade", "Clade", None), ("lineage", "Lineage", None), ("serotype", "Serotype", None),
           ("genotype", "Genotype", None), ("subtype", "Subtype", None), ("outbreak", "Outbreak", None),
           ("geographic_location", "Geographic Location", False), ("isolation_source", "Isolation Source", False),
           ("sra_accession", "SRA Accession", True), ("bioproject_accession", "BioProject", False),
           ("biosample_accession", "BioSample", False), ("sequencing_platform", "Sequencing Platform", False),
           ("assembly_method", "Assembly Method", False), ("has_sra", "SRA Available", False), ("has_fasta", "FASTA", False),
           ("genome_quality", "Genome Quality", False), ("genome_id", "Genome ID", False), ("bvbrc_genome_id", "BV-BRC ID", False),
           ("gc_content", "GC %", False), ("release_date", "Release Date", False), ("date_modified", "Modified", False)]
TAX_FACETS = [("clade", "Clade"), ("lineage", "Lineage"), ("serotype", "Serotype"), ("genotype", "Genotype"), ("subtype", "Subtype")]
GN_FACETS = [("genome_status", "Genome Status"), ("host_name", "Host"), ("isolation_country", "Country"),
             ("collection_year", "Collection Year"), ("has_sra", "SRA"), ("has_fasta", "FASTA"), ("genome_quality", "Quality")]


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
    status = Counter(r.get("genome_status") for r in g)
    return {
        "total": len(g), "complete_pct": round(100 * status.get("Complete", 0) / n, 1),
        "avg_len": round(sum(lens) / len(lens)) if lens else 0, "last_checked": last_checked,
        "pending": pending_total(), "prov": {p["source"]: p["n"] for p in prov},
        "monthly": by_prefix("release_date", 7), "year": by_prefix("collection_date", 4),
        "country": top("isolation_country"), "host": top("host_name", 8), "status": top("genome_status", 6),
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
    app = FastAPI(title="GenomeCurator")

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

    def render(request, name, **ctx):
        resp = T.TemplateResponse(request, name, {"badge": pending_count(), "watch": ctx.get("watch") or request.cookies.get("gc_watch", ""), **ctx})
        if ctx.get("watch"):
            resp.set_cookie("gc_watch", str(ctx["watch"]), max_age=86400 * 365, samesite="lax")
        return resp

    @app.get("/")
    def index(request: Request, msg: str = ""):
        return render(request, "index.html", msg=msg, nav="cur",
                      watches=_rows("SELECT w.*, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n_genomes FROM watchlist w"),
                      changesets=_rows("SELECT c.*, w.name wname FROM changeset c JOIN watchlist w ON w.id=c.watch_id "
                                       "WHERE c.status!='DISCARDED' ORDER BY c.id DESC LIMIT 30"))

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
        except ValueError as e:
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

    @app.get("/api/lookup")
    def api_lookup(q: str = ""):
        """Check a typed virus name against ICTV, NCBI Virus, BV-BRC and Pathoplexus (read-only)."""
        from .lookup import lookup_virus
        try:
            return lookup_virus(q)
        except ValueError as e:
            return JSONResponse({"error": str(e)}, 400)

    @app.post("/watch")
    def watch_add(name: str = Form(...), taxon: str = Form("")):
        try:
            curate.add_watch(name.strip(), taxon.strip() or None)
            msg = f"added {name}"
        except Exception as e:
            msg = f"error: {e}"
        return RedirectResponse(f"/?msg={msg}", 303)

    @app.post("/watch/{name}/check")
    def check_now(name: str):
        def work(p):
            cs = curate.check(name, progress=p)
            if not cs:
                return "/?msg=변경 사항 없음"
            curate.dryrun(cs, progress=p)
            return f"/changeset/{cs}"
        return start_job(work)

    @app.get("/api/job/{jid}")
    def job(jid: str):
        return JOBS.get(jid) or JSONResponse({"error": "unknown job", "done": True}, 404)

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
        agree = curate.source_agreement(cs["watch_id"])
        for i in _rows("SELECT * FROM changeset_item WHERE changeset_id=? ORDER BY kind, accession", cid):
            r = {k: v for k, v in json.loads(i["new_record_json"] or "{}").items()
                 if v is not None and k not in ("taxon_lineage_ids", "taxon_lineage_names")}
            if isinstance(r.get("submitter"), str):
                r["submitter"] = r["submitter"][:120]
            ag = agree.get(i["accession"], {})
            r.update(id=i["id"], kind=i["kind"], decision=i["decision"], accession=i["accession"],
                     has_sra="Yes" if r.get("sra_accession") else "No",
                     integrity="불일치" if ag.get("conflict") else ("표기차이" if ag.get("cosmetic") else "일치"),
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
        return start_job(work)

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
        return start_job(work)

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
        for n, r in enumerate(rows):
            r.pop("taxon_lineage_ids", None)
            r.pop("taxon_lineage_names", None)
            r["id"] = n
            r["genbank"] = f"{r['accession']}.{r['accession_version']}" if r.get("accession_version") else r["accession"]
            r["_prov"] = prov.get(r["genome_id"], {})
            r["has_sra"] = "Yes" if r.get("sra_accession") else "No"
            r["has_fasta"] = "Yes" if r.get("md5") else "No"
        ppx = _ppx_slug(watch, rows)
        data = {"rows": rows, "name": "genomes", "watch": cur and cur["name"], "lineage": lineage, "ppx": ppx,
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

    @app.get("/integrity")
    def integrity(request: Request, watch: int = 0, field: str = "", kind: str = "conflict", files: int = 0):
        from . import reconcile, verify
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
                      checks=verify.run(wid or None, check_files=bool(files)), files=files)

    @app.get("/dashboard")
    def dashboard(request: Request, watch: int = 0):
        watches = _rows("SELECT w.id, w.name, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n FROM watchlist w")
        watch = current_watch(request, watch, watches)
        g = _rows("SELECT * FROM genome WHERE watch_id=?", watch)
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
        return start_job(work)

    @app.get("/public")
    def public_view():
        """The same read-only page as docs/index.html, built live from the local database."""
        from .static_site import build_html
        try:
            return HTMLResponse(build_html("https://github.com/CHOICHI-HWAN/GenomeCurator"))
        except ValueError as e:
            return HTMLResponse(f"<p>{e}</p>", status_code=404)

    @app.get("/history")
    def history(request: Request, msg: str = ""):
        return render(request, "history.html", msg=msg, nav="his",
                      log=_rows("SELECT * FROM audit_log ORDER BY id DESC LIMIT 200"))

    if scheduler:
        from apscheduler.schedulers.background import BackgroundScheduler

        def sched_job():
            for w in _rows("SELECT name FROM watchlist WHERE enabled=1"):
                try:
                    cs = curate.check(w["name"])
                    if cs:
                        curate.dryrun(cs)
                except BaseException as e:
                    print("scheduled check failed:", w["name"], e)

        sched = BackgroundScheduler()
        sched.add_job(sched_job, "interval", hours=config.CHECK_INTERVAL_HOURS, id="check")
        sched.start()
        app.router.on_shutdown.append(lambda: sched.shutdown(wait=False))
    return app
