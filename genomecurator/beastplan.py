"""A whole dating analysis as one plan: temporal signal, then the clock x tree-prior comparison, then the dating itself.

Why a plan and not three separate runs: each step has to be judged on the same thing. Before this, the temporal signal and the
model comparison were run on a fresh 60-genome subsample while the dating used its own 166, and the substitution model was dropped
on the way into BETS so those steps silently ran HKY+gamma whatever the dating used. A conclusion about the clock then came from
data and a model that are not the ones the dates come from. A plan fixes the genomes, the alignment and the substitution model
once, writes them into `plan.json`, and every step reads them from there.

The substitution model is not searched here: ModelFinder chooses it (`beastprep.modelfinder_choice` reads the choice out of an
IQ-TREE report) and the plan takes it. What the plan compares, by marginal likelihood, is the clock and the tree prior - every
combination of them, which is `beastprep.CLOCKS` x `beastprep.PRIORS`.

The steps run as cluster jobs and the plan advances when a step's jobs are all in; `advance()` is called from the cluster's
background thread and is safe to call again at any time. A plan that finds no temporal signal stops there and says so, because a
tMRCA from data without a clock signal is a reading of the prior.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from . import beastprep, log

_log = log.get("beastplan")
STAGES = ("bets", "compare", "dating")
SIGNAL_MIN = 3.0                                                 # log Bayes factor for the dates: below this the plan stops (Kass & Raftery: 3 is "positive")


def plan_root(clock_run: dict) -> Path:
    return Path(clock_run["run_dir"]) / "beast" / "plan"


def plan_dir(clock_run: dict, name: str) -> Path:
    return plan_root(clock_run) / beastprep._safe(name)


def read(clock_run: dict, name: str) -> dict | None:
    f = plan_dir(clock_run, name) / "plan.json"
    try:
        return json.loads(f.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _write(clock_run: dict, p: dict) -> dict:
    d = plan_dir(clock_run, p["name"])
    d.mkdir(parents=True, exist_ok=True)
    (d / "plan.json").write_text(json.dumps(p, ensure_ascii=False, indent=1), encoding="utf-8")
    return p


def list_plans(clock_run: dict) -> list[dict]:
    r = plan_root(clock_run)
    return [json.loads(f.read_text(encoding="utf-8")) for f in sorted(r.glob("*/plan.json"))] if r.is_dir() else []


def start(clock_run: dict, name: str, tip_names: list[str], modelfinder: str | Path, dating: dict | None = None,
          bets: dict | None = None, chains: int = 2) -> dict:
    """Begin a plan. `tip_names` are the genomes every step uses, `modelfinder` the IQ-TREE report whose chosen model every step
    uses, `dating` the settings of the final analysis (chain length, priors, groups), `chains` how many full-length chains it runs."""
    from . import cluster
    if read(clock_run, name):
        raise ValueError(f"{name} 계획이 이미 있습니다.")
    if len(tip_names) < 10:
        raise ValueError("유전체가 10개 이상이어야 합니다.")
    model = beastprep.modelfinder_choice(modelfinder)
    tips = sorted(set(tip_names))                                      # one fixed, reproducible order for every step
    ds = dict(dating or {}, tip_names=tips, subst_model=model["subst_model"], site_rates=model["site_rates"])
    bs = dict(bets or {}, tip_names=tips, subst_model=model["subst_model"], site_rates=model["site_rates"], compare=False)
    p = {"name": name, "created": datetime.now().isoformat(timespec="seconds"), "run_dir": str(clock_run["run_dir"]),
         "clock_run_id": clock_run.get("id"), "n_tips": len(tips), "tip_names": tips, "model": model,
         "bets_settings": bs, "dating_settings": ds, "chains": int(chains), "stage": "bets", "groups": {}, "results": {},
         "grid": [f"{c}-{t}" for c in beastprep.CLOCKS for t in beastprep.PRIORS], "signal_min": SIGNAL_MIN}
    bs["name"] = f"plan-{beastprep._safe(name)}-signal"
    g = cluster.submit_bets(clock_run, bs)
    p["groups"]["bets"] = g["group"]
    p["log"] = [f"{p['created']}  시작: 유전체 {len(tip_names)}개, 치환 모형 {model['model']} ({model['note']}), 시계 신호 확인 {len(g['jobs'])}개 작업"]
    return _write(clock_run, p)


def _note(p: dict, msg: str) -> None:
    p.setdefault("log", []).append(f"{datetime.now().isoformat(timespec='seconds')}  {msg}")
    _log.info("plan %s: %s", p.get("name"), msg)


def advance(clock_run: dict, name: str) -> dict:
    """Move the plan on if its current step has finished. Safe to call again; does nothing when the step is still running."""
    from . import cluster
    p = read(clock_run, name)
    if not p or p["stage"] in ("done", "stopped"):
        return p or {}
    groups = {g["id"]: g for g in cluster.bets_groups()}
    splits = {g["id"]: g for g in cluster.split_groups()}

    if p["stage"] == "bets":
        g = groups.get(p["groups"].get("bets"))
        if not g or not g.get("joined"):
            return p
        res = g.get("summary") or {}
        bf = next((res.get("log_bf_" + k) for k in ("gss", "ss", "ps") if res.get("log_bf_" + k) is not None), None)
        p["results"]["bets"] = res
        if bf is None:
            p["stage"] = "stopped"
            p["stopped_because"] = "시계 신호의 베이즈 인수를 추정하지 못했습니다 (한계우도가 나오지 않음)."
            _note(p, p["stopped_because"])
            return _write(clock_run, p)
        if bf < p["signal_min"]:
            p["stage"] = "stopped"
            p["stopped_because"] = (f"수집일의 로그 베이즈 인수가 {bf:.1f}로 기준 {p['signal_min']}에 못 미칩니다. "
                                    "시계 신호가 없는 자료의 tMRCA는 사전분포를 되읽는 것에 가까우므로 여기서 멈춥니다.")
            _note(p, p["stopped_because"])
            return _write(clock_run, p)
        cs = dict(p["bets_settings"], compare=True, name=f"plan-{beastprep._safe(name)}-models")
        ng = cluster.submit_bets(clock_run, cs)
        p["groups"]["compare"] = ng["group"]
        p["compare_settings"] = cs
        p["stage"] = "compare"
        _note(p, f"시계 신호 확인: 로그 베이즈 인수 {bf:.1f} ≥ {p['signal_min']}. 모델 비교 {len(p['grid'])}조합, {len(ng['jobs'])}개 작업 시작")
        return _write(clock_run, p)

    if p["stage"] == "compare":
        g = groups.get(p["groups"].get("compare"))
        if not g or not g.get("joined"):
            return p
        res = g.get("summary") or {}
        ranking = (res.get("models") or {}).get("ranking") or []
        p["results"]["compare"] = res
        if not ranking:
            p["stage"] = "stopped"
            p["stopped_because"] = "모델 비교에서 한계우도를 얻지 못했습니다."
            _note(p, p["stopped_because"])
            return _write(clock_run, p)
        best = ranking[0]["model"]
        clock_m, _, tree_m = best.partition("-")
        p["chosen"] = {"clock": clock_m, "tree_prior": tree_m, "log_ml": ranking[0]["ml"],
                       "runner_up": ranking[1]["model"] if len(ranking) > 1 else None,
                       "log_bf_over_runner_up": -ranking[1]["log_bf_vs_best"] if len(ranking) > 1 else None}
        ds = dict(p["dating_settings"], clock=clock_m, tree_prior=tree_m, tag=f"plan-{beastprep._safe(name)}-dating")
        beastprep.prepare(clock_run, ds)                                   # the dating package, under the model the comparison chose
        d = beastprep.folder(clock_run, ds["tag"])
        ng = cluster.submit_beast_split(d, p["chains"], f"plan {name}: dating under {best}", mode="replicate")
        p["groups"]["dating"] = ng["group"]
        p["dating_settings"] = ds
        p["dating_out"] = ng["out"]
        p["stage"] = "dating"
        gap = p["chosen"]["log_bf_over_runner_up"]
        _note(p, f"모델 비교 1위 {best}" + (f" (2위 {p['chosen']['runner_up']}와 로그 베이즈 인수 {gap:.1f})" if gap is not None else "")
              + f". 본 분석 {p['chains']}개 사슬 시작")
        return _write(clock_run, p)

    if p["stage"] == "dating":
        g = splits.get(p["groups"].get("dating"))
        if not g or not g.get("joined"):
            return p
        out = Path(p["dating_out"])
        p["results"]["dating"] = {"folder": str(out), "agreement": g.get("agreement")}
        p["stage"] = "done"
        ok = (g.get("agreement") or {}).get("ok")
        _note(p, f"본 분석 완료, 사슬 {'일치' if ok else '불일치 — 확인 필요'}. 결과 {out.name}")
        return _write(clock_run, p)
    return p


def unfinished() -> list[dict]:
    """Every plan on disk that has not finished, newest first. Found by looking for plan.json, so a plan survives a restart."""
    from . import config
    out = []
    for f in sorted(config.DATA_DIR.glob("nextstrain/runs/*/*/beast/plan/*/plan.json"), key=lambda x: x.stat().st_mtime, reverse=True):
        try:
            p = json.loads(f.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if p.get("stage") not in ("done", "stopped"):
            out.append(p)
    return out


def advance_all() -> list[dict]:
    """Move on every plan whose step has finished. Called from the cluster's background thread; never raises."""
    from . import nextstrain
    out = []
    for p in unfinished():
        try:
            cr = nextstrain.get_run(p["clock_run_id"]) if p.get("clock_run_id") else None
            if not cr:
                _log.warning("plan %s: its clock run %s is gone", p.get("name"), p.get("clock_run_id"))
                continue
            out.append(advance(cr, p["name"]))
        except Exception as ex:
            _log.warning("plan %s could not be advanced: %s", p.get("name"), ex)
    return out


def status(clock_run: dict, name: str) -> dict:
    """What the plan looks like now, with the jobs of its current step."""
    from . import cluster
    p = read(clock_run, name)
    if not p:
        return {}
    gid = p["groups"].get(p["stage"]) if p["stage"] in STAGES else None
    groups = {g["id"]: g for g in (cluster.bets_groups() + cluster.split_groups())}
    g = groups.get(gid) if gid else None
    return {**p, "current_group": g, "done_jobs": (g or {}).get("done"), "total_jobs": (g or {}).get("expected") or (g or {}).get("parts")}
