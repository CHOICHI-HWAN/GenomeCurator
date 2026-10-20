"""All the analyses of a dataset in one go: QC -> Squirrel -> mutation profile -> variant calls (snippy + lofreq) -> Nextstrain, one after the other as one background job.

Each step is the same function the single-step buttons call, so nothing is computed differently here: a step whose input is unchanged (same members, content, settings and tool version) is reused
instead of run again, and a step that fails does not stop the others (the variant calls are skipped when the mutation profile they belong to failed). The result says, per step, whether it was run,
reused, failed or skipped, and why. Nextstrain is last because it takes longest; the quick ones are finished and usable while it runs."""
from . import clock, datasets, mutprofile, nextstrain, squirrel, support, variantflow

STEPS = ("qc", "squirrel", "clock", "mutprofile", "variants", "nextstrain")
LABEL = {"qc": "QC", "squirrel": "Squirrel", "clock": "분자시계·tMRCA", "mutprofile": "변이 프로파일", "variants": "변이 호출(snippy+lofreq)", "nextstrain": "Nextstrain 트리"}
WEIGHT = {"qc": 5, "squirrel": 25, "clock": 4, "mutprofile": 12, "variants": 8, "nextstrain": 46}


def plan(ds: dict, steps=None) -> tuple[list[str], list[str]]:
    """(steps to run in order, notes). The analyses are built for monkeypox only (support.py): for another virus only QC is run."""
    want = [s for s in STEPS if (steps is None or s in steps)]
    notes = []
    from .db import connect
    con = connect()
    try:
        w = con.execute("SELECT name, taxon FROM watchlist WHERE id=?", (ds["watch_id"],)).fetchone()
    finally:
        con.close()
    level = {s["area"]: s["level"] for s in support.for_virus(w["name"], w["taxon"])} if w else {}
    if level.get("analysis") == "none":
        dropped = [s for s in want if s != "qc"]
        if dropped:
            notes.append("이 바이러스는 분석(Squirrel·변이·Nextstrain)이 아직 준비되지 않아 QC만 실행합니다.")
        want = [s for s in want if s == "qc"]
    if "variants" in want and "mutprofile" not in want:
        notes.append("변이 호출은 변이 프로파일 결과가 필요한데 이번에 고르지 않아, 이미 있는 가장 최근 변이 프로파일 결과에 돌립니다.")
    return want, notes


def run_all(ds: dict, steps=None, force: bool = False, cores: int | None = None, progress=None) -> dict:
    p = progress or (lambda *a: None)
    want, notes = plan(ds, steps)
    total = sum(WEIGHT[s] for s in want) or 1
    done_w = 0
    out: dict = {"steps": {}, "notes": notes}
    mp_run = None

    def sub(step):                                      # progress of one step -> the share of the whole it owns
        base = done_w
        return lambda pct, m="": p((base + WEIGHT[step] * float(pct) / 100) * 100 / total, f"[{LABEL[step]}] {m}")

    for s in want:
        sp = sub(s)
        try:
            if s == "qc":
                res = datasets.run_qc(ds, sp, force)
                out["steps"][s] = {"status": "done", "detail": f"공개 {((res.get('public') or {}).get('sequences') or 0)}건"}
            elif s == "squirrel":
                r = squirrel.run(ds, cores, force, "전체 분석", sp)
                out["steps"][s] = {"status": "reused" if r.get("reused") else "done", "run": r["id"]}
            elif s == "clock":
                if not clock.find_source(ds, "squirrel") and not clock.find_source(ds, "nextstrain"):
                    out["steps"][s] = {"status": "skipped", "detail": "트리가 있는 Squirrel 또는 Nextstrain 결과가 없습니다"}
                else:
                    r = clock.run(ds, {"source": "squirrel" if clock.find_source(ds, "squirrel") else "nextstrain"}, force, "전체 분석", sp)
                    out["steps"][s] = {"status": "reused" if r.get("reused") else "done", "run": r["id"]}
            elif s == "mutprofile":
                r = mutprofile.run(ds, cores, force, "전체 분석", sp)
                mp_run = nextstrain.get_run(r["id"])
                out["steps"][s] = {"status": "reused" if r.get("reused") else "done", "run": r["id"]}
            elif s == "variants":
                r = mp_run or next((x for x in nextstrain.list_runs(ds["id"]) if x["tool"] == "mutprofile" and x["status"] == "done"), None)
                if not r:
                    out["steps"][s] = {"status": "skipped", "detail": "변이 프로파일 결과가 없습니다"}
                elif mp_run is None and "mutprofile" in want:
                    out["steps"][s] = {"status": "skipped", "detail": "변이 프로파일이 실패했습니다"}
                elif variantflow.available(r) and not force:
                    out["steps"][s] = {"status": "reused", "run": r["id"], "detail": "이 변이 프로파일 결과에는 이미 호출 결과가 있습니다"}
                else:
                    x = variantflow.run_recorded(r, cores, sp)
                    out["steps"][s] = {"status": "done", "run": x.get("run"), "detail": f"캐시에서 {x.get('from_cache', 0)}건, 새로 {x.get('genomes', 0) - x.get('from_cache', 0)}건"}
            elif s == "nextstrain":
                r = nextstrain.run(ds, "ambient", cores, force, "전체 분석", sp)
                out["steps"][s] = {"status": "reused" if r.get("reused") else "done", "run": r["id"]}
        except Exception as e:                           # one failing step must not stop the others
            out["steps"][s] = {"status": "error", "detail": str(e)[:300]}
            if s == "mutprofile":
                mp_run = None
        done_w += WEIGHT[s]
    ok = sum(1 for v in out["steps"].values() if v["status"] in ("done", "reused"))
    p(100, f"{ok}/{len(want)}단계 완료" + (f", 실패 {sum(1 for v in out['steps'].values() if v['status'] == 'error')}" if any(v["status"] == "error" for v in out["steps"].values()) else ""))
    return out
