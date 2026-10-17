"""Background jobs: every long task (adding a virus, check, apply, FASTA download, sequence verification, schema check, publish,
and the scheduler's own runs) runs on the server as a job with a durable record.

* live state lives in memory (fast to poll); every change is also written to the `job` table, so a job can still be looked at after a page
  reload, after the live entry is gone, and after a server restart
* a job that was running when the server stopped is marked *interrupted* at the next start (it will not silently look "running" forever)
* each job keeps its own timestamped log, the error with a traceback tail, and optionally a result (e.g. what a virus lookup found)
The web UI (`/jobs`, the progress modal, the badge and toast in every page) and the scheduler both go through here."""
import json
import threading
import time
import traceback
import uuid
from datetime import datetime
from pathlib import Path

from . import log
from .db import audit_now, connect

_log = log.get("jobs")
live: dict[str, dict] = {}
KEEP = 300                     # finished jobs kept in the table
LOG_MAX = 800
SAVE_EVERY = 1.0               # seconds between progress writes (the first and the last are always written)
KINDS = {"check": "확인", "apply": "적용", "fasta": "서열 받기", "seqcheck": "서열 검증", "schema": "스키마 점검", "watch_add": "바이러스 추가",
         "lookup": "바이러스 조회", "deploy": "배포", "scheduled": "예약 실행", "job": "작업"}


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _save(j: dict) -> None:
    try:
        con = connect()
        try:
            con.execute("INSERT OR REPLACE INTO job(id,kind,label,status,pct,msg,started,finished,error,redirect,log_json,result_json) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                        (j["id"], j["kind"], j["label"], j["status"], j["pct"], j["msg"], j["started"], j.get("finished"), j.get("error"), j.get("redirect"),
                         json.dumps(j["log"], ensure_ascii=False), json.dumps(j["result"], ensure_ascii=False) if j.get("result") is not None else None))
            con.commit()
        finally:
            con.close()
        j["_saved"] = time.monotonic()
    except Exception:                   # a bookkeeping failure must never break the job itself
        _log.warning("could not save job %s", j.get("id"), exc_info=True)


def _prune() -> None:
    try:
        con = connect()
        try:
            con.execute("DELETE FROM job WHERE status!='running' AND id NOT IN (SELECT id FROM job ORDER BY started DESC LIMIT ?)", (KEEP,))
            con.commit()
        finally:
            con.close()
    except Exception:
        _log.debug("prune failed", exc_info=True)


def _public(j: dict) -> dict:
    out = {k: j.get(k) for k in ("id", "kind", "label", "status", "pct", "msg", "started", "finished", "error", "redirect", "log", "result")}
    out["done"] = j["status"] != "running"
    out["kind_label"] = KINDS.get(j["kind"], j["kind"])
    return out


def run(jid: str, label: str, kind: str, fn) -> dict:
    """Run fn(progress) in the CURRENT thread as job `jid`. fn returns a redirect URL (str) or {"redirect":..., "result":...}."""
    j = {"id": jid, "kind": kind, "label": label, "status": "running", "pct": 0, "msg": "시작 중…", "started": _now(), "finished": None,
         "error": None, "redirect": "/", "log": [], "result": None, "_last": None, "_saved": 0.0}
    live[jid] = j

    def addlog(msg: str, level: str = "info"):
        if len(j["log"]) < LOG_MAX:
            j["log"].append({"t": time.strftime("%H:%M:%S"), "msg": msg, "level": level})

    def prog(pct, msg):
        j["pct"], j["msg"] = pct, msg
        if msg != j["_last"]:
            j["_last"] = msg
            addlog(f"{msg}" + (f"  ({round(pct)}%)" if 1 < pct < 100 else ""))
            _log.info("job %s progress %s%%: %s", label, round(pct), msg)
            if time.monotonic() - j["_saved"] >= SAVE_EVERY:
                _save(j)

    log.set_context(job=jid)
    t0 = time.monotonic()
    _log.info("job %s started (%s)", label, kind)
    addlog(f"작업 시작: {label}")
    _save(j)
    try:
        res = fn(prog)
        if isinstance(res, dict):
            j["redirect"], j["result"] = res.get("redirect"), res.get("result")
        elif res is not None:
            j["redirect"] = res
        j["pct"], j["status"] = 100, "done"
        addlog("완료")
        _log.info("job %s finished in %.1fs", label, time.monotonic() - t0)
    except BaseException as e:
        msg = str(e) or e.__class__.__name__
        j["error"], j["status"] = msg, "error"
        tb = traceback.extract_tb(e.__traceback__)
        where = f"{Path(tb[-1].filename).name}:{tb[-1].lineno} in {tb[-1].name}" if tb else ""
        addlog(f"오류: {type(e).__name__}: {msg}" + (f"  [{where}]" if where else ""), "error")
        for fr in tb[-4:]:      # the tail of the call chain, right here in the log the person reads
            addlog(f"   {Path(fr.filename).name}:{fr.lineno} {fr.name}", "error")
        addlog("서버 로그에 전체 traceback이 있습니다 (이력 → 서버 로그).", "error")
        _log.error("job %s failed after %.1fs: %s: %s", label, time.monotonic() - t0, type(e).__name__, msg,
                   exc_info=not isinstance(e, SystemExit) and type(e).__name__ != "SourceUnavailable")
        audit_now("job_failed", f"{label} ({jid}): {type(e).__name__}: {msg}")
    j["finished"] = _now()
    _save(j)
    _prune()
    out = _public(j)
    # keep finished jobs in memory only briefly: the table is the record
    t = threading.Timer(120, lambda: live.pop(jid, None))
    t.daemon = True            # must not hold the process open at exit
    t.start()
    return out


def start(fn, label: str, kind: str = "job") -> str:
    """Run fn on a background thread; returns the job id at once."""
    jid = uuid.uuid4().hex[:8]
    live[jid] = {"id": jid, "kind": kind, "label": label, "status": "running", "pct": 0, "msg": "시작 중…", "started": _now(), "finished": None,
                 "error": None, "redirect": "/", "log": [], "result": None, "_last": None, "_saved": 0.0}   # visible to pollers before the thread begins
    threading.Thread(target=run, args=(jid, label, kind, fn), daemon=True, name=f"job-{jid}").start()
    return jid


def run_sync(label: str, kind: str, fn) -> dict:
    """For the scheduler: the same record and log, but in the calling thread."""
    return run(uuid.uuid4().hex[:8], label, kind, fn)


def _row(r) -> dict:
    d = dict(r)
    d["log"] = json.loads(d.pop("log_json") or "[]")
    rj = d.pop("result_json", None)
    d["result"] = json.loads(rj) if rj else None
    return _public(d)


def get(jid: str) -> dict | None:
    if jid in live:
        return _public(live[jid])
    con = connect()
    try:
        r = con.execute("SELECT * FROM job WHERE id=?", (jid,)).fetchone()
        return _row(r) if r else None
    finally:
        con.close()


def list_jobs(limit: int = 100, status: str = "", kind: str = "") -> list[dict]:
    """Newest first, without the logs (they can be long); running jobs come from memory so their progress is current."""
    q, a = "SELECT id,kind,label,status,pct,msg,started,finished,error,redirect FROM job WHERE 1=1", []
    if status:
        q += " AND status=?"
        a.append(status)
    if kind:
        q += " AND kind=?"
        a.append(kind)
    q += " ORDER BY started DESC, rowid DESC LIMIT ?"
    a.append(limit)
    con = connect()
    try:
        rows = [dict(r) for r in con.execute(q, a)]
    finally:
        con.close()
    seen = {r["id"] for r in rows}
    for jid, j in list(live.items()):
        if j["status"] == "running" and (not status or status == "running") and (not kind or kind == j["kind"]):
            if jid not in seen:
                rows.append({k: j.get(k) for k in ("id", "kind", "label", "status", "pct", "msg", "started", "finished", "error", "redirect")})
        for r in rows:
            if r["id"] == jid:
                r.update({k: j.get(k) for k in ("status", "pct", "msg", "finished", "error")})
    rows.sort(key=lambda r: r["started"] or "", reverse=True)
    for r in rows:
        r["kind_label"] = KINDS.get(r["kind"], r["kind"])
        r["done"] = r["status"] != "running"
    return rows[:limit]


def summary(since: str = "", fail_since: str = "") -> dict:
    """What every page needs for its badge and toast: how many are running, and what finished (or failed) since `since`."""
    rows = list_jobs(60)
    running = [{"id": r["id"], "label": r["label"], "pct": r["pct"], "msg": r["msg"]} for r in rows if r["status"] == "running"]
    fin = [r for r in rows if r["status"] != "running" and r.get("finished") and (not since or r["finished"] > since)]
    unseen = sum(1 for r in rows if r["status"] in ("error", "interrupted") and r.get("finished") and (not fail_since or r["finished"] > fail_since))
    return {"running": running, "failed": sum(1 for r in fin if r["status"] in ("error", "interrupted")), "failed_unseen": unseen,
            "finished": [{"id": r["id"], "label": r["label"], "status": r["status"], "error": r["error"], "finished": r["finished"], "redirect": r.get("redirect")} for r in fin[:8]],
            "now": _now()}


def mark_interrupted() -> int:
    """At server start: whatever the table still calls running died with the previous process."""
    con = connect()
    try:
        n = con.execute("UPDATE job SET status='interrupted', finished=?, error=? WHERE status='running'",
                        (_now(), "서버가 다시 시작되어 중단되었습니다. 필요하면 같은 작업을 다시 실행하세요.")).rowcount
        con.commit()
    finally:
        con.close()
    if n:
        _log.warning("%d job(s) were running when the server stopped; marked as interrupted", n)
    return n
