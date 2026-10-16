import json
from pathlib import Path

from fastapi import FastAPI, Form, Request
from fastapi.responses import RedirectResponse, Response
from fastapi.templating import Jinja2Templates

from . import config, curate
from .db import connect

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


def create_app(scheduler: bool = True) -> FastAPI:
    app = FastAPI(title="GenomeCurator")

    def pending_count():
        return _rows("SELECT COUNT(*) n FROM changeset WHERE status IN ('READY','BLOCKED')")[0]["n"]

    def render(request, name, **ctx):
        return T.TemplateResponse(request, name, {"badge": pending_count(), **ctx})

    @app.get("/")
    def index(request: Request, msg: str = ""):
        return render(request, "index.html", msg=msg,
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
        try:
            cs = curate.check(name)
            if cs:
                curate.dryrun(cs)
                return RedirectResponse(f"/changeset/{cs}", 303)
            return RedirectResponse("/?msg=no changes", 303)
        except BaseException as e:
            return RedirectResponse(f"/?msg=check failed: {e}", 303)

    @app.get("/changeset/{cid}")
    def changeset(request: Request, cid: int):
        cs = _rows("SELECT c.*, w.name wname FROM changeset c JOIN watchlist w ON w.id=c.watch_id WHERE c.id=?", cid)[0]
        items = _rows("SELECT * FROM changeset_item WHERE changeset_id=? ORDER BY kind, accession", cid)
        for i in items:
            i["diff"] = json.loads(i["diff_json"] or "{}")
            i["prov"] = json.loads(i["provenance_json"] or "{}")
            i["rec"] = json.loads(i["new_record_json"] or "{}")
        return render(request, "changeset.html", cs=cs, items=items,
                      summary=json.loads(cs["summary_json"] or "{}"),
                      dry=json.loads(cs["dryrun_json"] or "{}"))

    @app.post("/changeset/{cid}/decide")
    async def decide(request: Request, cid: int):
        form = await request.form()
        ids = [int(x) for x in form.getlist("item")]
        action = form.get("action")
        if action in ("approved", "rejected", "hold"):
            curate.decide(ids, action)
        elif action == "approve_all":
            curate.decide([i["id"] for i in _rows(
                "SELECT id FROM changeset_item WHERE changeset_id=? AND decision='pending'", cid)], "approved")
        return RedirectResponse(f"/changeset/{cid}", 303)

    @app.post("/changeset/{cid}/approve_apply")
    def approve_apply(cid: int):
        curate.decide([i["id"] for i in _rows(
            "SELECT id FROM changeset_item WHERE changeset_id=? AND decision='pending'", cid)], "approved")
        return do_apply(cid)

    @app.get("/changeset/{cid}.tsv")
    def changeset_tsv(cid: int):
        rows = _rows("SELECT * FROM changeset_item WHERE changeset_id=? ORDER BY kind, accession", cid)
        out = []
        for r in rows:
            rec = json.loads(r["new_record_json"] or "{}")
            diff = json.loads(r["diff_json"] or "{}")
            out.append({"kind": r["kind"], "decision": r["decision"], "accession": r["accession"], **rec,
                        "changed_fields": ";".join(diff)})
        return _tsv(out, f"changeset_{cid}.tsv")

    @app.get("/genomes.tsv")
    def genomes_tsv():
        return _tsv(_rows("SELECT g.*, s.gc_content, s.md5, s.file_path FROM genome g "
                          "LEFT JOIN genome_sequence s ON s.accession=g.accession ORDER BY g.genome_id"), "genomes.tsv")

    @app.post("/changeset/{cid}/dryrun")
    def rerun(cid: int):
        curate.dryrun(cid)
        return RedirectResponse(f"/changeset/{cid}", 303)

    @app.post("/changeset/{cid}/apply")
    def do_apply(cid: int):
        try:
            res = curate.apply(cid)
            msg = f"applied {res['applied']} items"
        except BaseException as e:
            msg = f"apply failed: {e}"
        return RedirectResponse(f"/history?msg={msg}", 303)

    @app.get("/genomes")
    def genomes(request: Request):
        return render(request, "genomes.html", genomes=_rows(
            "SELECT g.*, s.md5, s.file_path FROM genome g LEFT JOIN genome_sequence s ON s.accession=g.accession ORDER BY g.genome_id"))

    @app.get("/history")
    def history(request: Request, msg: str = ""):
        return render(request, "history.html", msg=msg, log=_rows("SELECT * FROM audit_log ORDER BY id DESC LIMIT 200"))

    if scheduler:
        from apscheduler.schedulers.background import BackgroundScheduler

        def job():
            for w in _rows("SELECT name FROM watchlist WHERE enabled=1"):
                try:
                    cs = curate.check(w["name"])
                    if cs:
                        curate.dryrun(cs)
                except BaseException as e:
                    print("scheduled check failed:", w["name"], e)

        sched = BackgroundScheduler()
        sched.add_job(job, "interval", hours=config.CHECK_INTERVAL_HOURS, id="check")
        sched.start()
        app.router.on_shutdown.append(lambda: sched.shutdown(wait=False))
    return app
