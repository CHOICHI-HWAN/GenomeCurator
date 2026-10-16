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
GN_COLS = [("genome_name", "Genome Name", True), ("strain", "Strain", True), ("genbank", "GenBank Accessions", True),
           ("genome_length", "Size", True), ("genome_status", "Genome Status", True), ("collection_date", "Collection Date", True),
           ("isolation_country", "Isolation Country", True), ("host_name", "Host Name", True),
           ("geographic_location", "Geographic Location", False), ("isolation_source", "Isolation Source", False),
           ("genome_quality", "Genome Quality", False), ("genome_id", "Genome ID", False), ("bvbrc_genome_id", "BV-BRC ID", False),
           ("gc_content", "GC %", False), ("release_date", "Release Date", False), ("date_modified", "Modified", False)]
GN_FACETS = [("genome_status", "Genome Status"), ("host_name", "Host"), ("isolation_country", "Country"),
             ("collection_year", "Collection Year"), ("genome_quality", "Quality")]


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
              "geographic_location", "genotype", "subtype", "genome_quality", "publication"]
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
        "hist": hist, "fill": {"labels": [x["f"] for x in fill], "values": [x["pct"] for x in fill]},
    }


def create_app(scheduler: bool = True) -> FastAPI:
    app = FastAPI(title="GenomeCurator")

    def pending_count():
        return pending_total()

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
            r = {k: v for k, v in json.loads(i["new_record_json"] or "{}").items()
                 if v is not None and k not in ("taxon_lineage_ids", "taxon_lineage_names")}
            if isinstance(r.get("submitter"), str):
                r["submitter"] = r["submitter"][:120]
            r.update(id=i["id"], kind=i["kind"], decision=i["decision"], accession=i["accession"])
            if i["kind"] == "changed":
                r["_diff"] = json.loads(i["diff_json"] or "{}")
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
        action = form.get("decision")
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
    def genomes(request: Request, watch: int = 0):
        watches = _rows("SELECT w.id, w.name, w.taxon_id, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n FROM watchlist w")
        if not watch and watches:
            watch = next((w["id"] for w in watches if w["n"]), watches[0]["id"])
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
        data = {"rows": rows, "name": "genomes", "watch": cur and cur["name"], "lineage": lineage,
                "cols": [{"k": k, "label": l, "on": on} for k, l, on in GN_COLS],
                "facets": [{"k": k, "label": l} for k, l in GN_FACETS]}
        return render(request, "genomes.html", data=data, nav="gen", watches=watches, watch=watch)

    @app.post("/genomes/fasta")
    async def genomes_fasta(request: Request):
        form = await request.form()
        ids = [i for i in str(form.get("ids", "")).split(",") if i]
        out = []
        for r in _rows("SELECT file_path FROM genome_sequence WHERE genome_id IN (%s)" % ",".join("?" * len(ids)), *ids) if ids else []:
            try:
                out.append(Path(r["file_path"]).read_text())
            except OSError:
                pass
        return Response("".join(out), media_type="text/plain; charset=utf-8",
                        headers={"Content-Disposition": "attachment; filename=genomes.fna"})

    @app.get("/genomes.tsv")
    def genomes_tsv():
        return _tsv(_rows("SELECT g.*, s.gc_content, s.md5, s.file_path FROM genome g "
                          "LEFT JOIN genome_sequence s ON s.accession=g.accession ORDER BY g.genome_id"), "genomes.tsv")

    @app.get("/dashboard")
    def dashboard(request: Request, watch: int = 0):
        watches = _rows("SELECT w.id, w.name, (SELECT COUNT(*) FROM genome g WHERE g.watch_id=w.id) n FROM watchlist w")
        if not watch and watches:
            watch = next((w["id"] for w in watches if w["n"]), watches[0]["id"])
        g = _rows("SELECT * FROM genome WHERE watch_id=?", watch)
        prov = _rows("SELECT p.source, COUNT(*) n FROM field_provenance p JOIN genome g ON g.genome_id=p.genome_id "
                     "WHERE g.watch_id=? GROUP BY p.source", watch)
        last = _rows("SELECT last_checked FROM watchlist WHERE id=?", watch)
        return render(request, "dashboard.html", nav="dash", watches=watches, watch=watch,
                      data=dashboard_data(g, prov, last[0]["last_checked"] if last else None))

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
