import json
import threading
import uuid
from pathlib import Path

from fastapi import FastAPI, Form, Request
from fastapi.responses import JSONResponse, RedirectResponse, Response
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


CS_COLS = [("kind", "종류", True, True), ("decision", "결정", True, True), ("accession", "Accession", True, False),
           ("accession_version", "Ver", True, False), ("genome_name", "Genome Name", True, False),
           ("genome_status", "Status", True, False), ("genome_length", "Length", True, False),
           ("strain", "Strain", True, False), ("host_name", "Host", True, False),
           ("isolation_country", "Country", True, False), ("geographic_location", "Geographic Location", False, False),
           ("collection_date", "Collection Date", True, False), ("isolation_source", "Isolation Source", False, False),
           ("submitter", "Submitter", False, False), ("bvbrc_genome_id", "BV-BRC ID", False, False),
           ("release_date", "Release Date", False, False)]
CS_FACETS = [("kind", "종류"), ("decision", "결정"), ("genome_status", "Genome Status"), ("host_name", "Host"),
             ("isolation_country", "Country"), ("collection_year", "Collection Year")]
GN_COLS = [("genome_id", "Genome ID", True), ("accession", "Accession", True), ("genome_name", "Genome Name", True),
           ("strain", "Strain", True), ("genome_status", "Status", True), ("genome_length", "Length", True),
           ("host_name", "Host", True), ("isolation_country", "Country", True),
           ("geographic_location", "Geographic Location", False), ("collection_date", "Collection Date", True),
           ("isolation_source", "Isolation Source", False), ("genome_quality", "Quality", False),
           ("bvbrc_genome_id", "BV-BRC ID", False), ("md5", "MD5", False), ("date_modified", "Modified", False)]
GN_FACETS = [("genome_status", "Genome Status"), ("host_name", "Host"), ("isolation_country", "Country"),
             ("collection_year", "Collection Year"), ("genome_quality", "Quality")]


def create_app(scheduler: bool = True) -> FastAPI:
    app = FastAPI(title="GenomeCurator")

    def pending_count():
        return _rows("SELECT COUNT(*) n FROM changeset WHERE status IN ('READY','BLOCKED')")[0]["n"]

    def render(request, name, **ctx):
        return T.TemplateResponse(request, name, {"badge": pending_count(), **ctx})

    @app.get("/")
    def index(request: Request, msg: str = ""):
        return render(request, "index.html", msg=msg, nav="cur",
                      watches=_rows("SELECT w.*, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n_genomes FROM watchlist w"),
                      changesets=_rows("SELECT c.*, w.name wname FROM changeset c JOIN watchlist w ON w.id=c.watch_id "
                                       "WHERE c.status!='DISCARDED' ORDER BY c.id DESC LIMIT 30"))

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

    @app.get("/changeset/{cid}")
    def changeset(request: Request, cid: int):
        cs = _rows("SELECT c.*, w.name wname FROM changeset c JOIN watchlist w ON w.id=c.watch_id WHERE c.id=?", cid)[0]
        rows = []
        for i in _rows("SELECT * FROM changeset_item WHERE changeset_id=? ORDER BY kind, accession", cid):
            r = {k: v for k, v in json.loads(i["new_record_json"] or "{}").items() if v is not None}
            r.update(id=i["id"], kind=i["kind"], decision=i["decision"], accession=i["accession"],
                     _diff=json.loads(i["diff_json"] or "{}"))
            rows.append(r)
        data = {"rows": rows, "selectable": True, "name": f"changeset_{cid}",
                "cols": [{"k": k, "label": l, "on": on, "cls": c} for k, l, on, c in CS_COLS],
                "facets": [{"k": k, "label": l} for k, l in CS_FACETS]}
        return render(request, "changeset.html", cs=cs, data=data, nav="cur",
                      summary=json.loads(cs["summary_json"] or "{}"),
                      dry=json.loads(cs["dryrun_json"] or "{}"))

    @app.post("/changeset/{cid}/decide")
    async def decide(request: Request, cid: int):
        form = await request.form()
        ids = [int(x) for x in form.getlist("item")]
        action = form.get("action")
        if action in ("approved", "rejected", "hold"):
            curate.decide(ids, action)
        return RedirectResponse(f"/changeset/{cid}", 303)

    def _apply_job(cid: int):
        def work(p):
            res = curate.apply(cid, progress=p)
            return f"/history?msg={res['applied']}건 적용 완료"
        return start_job(work)

    @app.post("/changeset/{cid}/approve_apply")
    async def approve_apply(request: Request, cid: int):
        form = await request.form()
        curate.decide([int(x) for x in form.getlist("item")], "approved")
        return _apply_job(cid)

    @app.post("/changeset/{cid}/apply")
    def do_apply(cid: int):
        return _apply_job(cid)

    @app.post("/changeset/{cid}/dryrun")
    def rerun(cid: int):
        curate.dryrun(cid)
        return RedirectResponse(f"/changeset/{cid}", 303)

    @app.get("/genomes")
    def genomes(request: Request):
        rows = _rows("SELECT g.*, s.md5 FROM genome g LEFT JOIN genome_sequence s ON s.accession=g.accession ORDER BY g.genome_id")
        for n, r in enumerate(rows):
            r["id"] = n
            if r.get("accession_version"):
                r["accession"] = f"{r['accession']}.{r['accession_version']}"
            r["md5"] = (r.get("md5") or "")[:8]
        data = {"rows": rows, "selectable": False, "name": "genomes",
                "cols": [{"k": k, "label": l, "on": on} for k, l, on in GN_COLS],
                "facets": [{"k": k, "label": l} for k, l in GN_FACETS]}
        return render(request, "genomes.html", data=data, nav="gen")

    @app.get("/genomes.tsv")
    def genomes_tsv():
        return _tsv(_rows("SELECT g.*, s.gc_content, s.md5, s.file_path FROM genome g "
                          "LEFT JOIN genome_sequence s ON s.accession=g.accession ORDER BY g.genome_id"), "genomes.tsv")

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
