"""Molecular clock and tMRCA of a dataset: a root-to-tip regression (what TempEst does) and a time-scaled tree with confidence intervals (TreeTime), from the dataset's tree and the collection dates.

**Stage 1 - is there a clock, and how fast?** `treetime clock` regresses every tip's distance from the root on its collection date (TempEst's method; TempEst itself is a Java GUI that cannot be run on
a server, TreeTime does the same regression and also finds the best-fit root and the tips that do not fit). Result: the substitution rate, R^2 (how clock-like the tree is) and the root date
(where the line meets zero distance), plus the outlier tips.

**Stage 2 - when did the clades start?** `treetime` (branch lengths of the tree kept as they are) turns the tree into a time-scaled tree: a date with a 90% interval for every node. The dates of the root
and of the clades / lineages that are monophyletic in the tree (the MRCA of their genomes) are reported - the tMRCA.

The tree is the one the dataset already has: **Squirrel's** maximum-likelihood tree (default) or the **Nextstrain** tree (its divergence in substitutions, divided by the genome length). Dates are the curated
collection dates: a full date is used as it is, a month as the range of that month (a year alone is not used unless asked for), a missing date drops the genome.

What this cannot tell: the root date is an extrapolation from the sampled period, and a weak clock (low R^2, a short period, rate dominated by APOBEC3 changes) gives a very wide interval whatever the point
estimate says - `warnings()` states this next to the numbers. For a proper Bayesian estimate see `beastprep.py`."""
import csv
import hashlib
import json
import math
import re
import subprocess
from collections import Counter
from datetime import date, datetime
from pathlib import Path

from . import aligned, datasets, lineages, log, nextstrain, runs
from .db import connect

_log = log.get("clock")
DEFAULTS = {"source": "squirrel", "min_precision": "month", "clock_filter": 3.0, "reroot": "least-squares", "covariation": True, "coalescent": "", "min_group": 5}   # covariation: stage 2 only
SOURCES = ("squirrel", "nextstrain")
PRECISION = ("day", "month", "year")
REROOTS = ("least-squares", "min_dev", "oldest", "keep")


def env_prefix() -> Path | None:
    pre = nextstrain.env_prefix()
    return pre if pre and (pre / "bin" / "treetime").exists() else None


def status() -> dict:
    pre = env_prefix()
    ver = ""
    if pre:
        try:
            ver = subprocess.run([str(pre / "bin" / "treetime"), "--version"], capture_output=True, text=True, timeout=60).stdout.strip().split()[-1]
        except Exception:
            ver = ""
    return {"ok": bool(pre), "version": ver, "how": "mamba env create -f nextstrain/environment.yml   (TreeTime is part of the Nextstrain environment)"}


def parse_settings(raw: dict | None) -> dict:
    s = {**DEFAULTS, **(raw or {})}
    out = {"source": str(s["source"]), "min_precision": str(s["min_precision"]), "reroot": str(s["reroot"]), "covariation": bool(s["covariation"]), "coalescent": str(s["coalescent"] or "").strip()}
    if out["source"] not in SOURCES:
        raise ValueError(f"트리 출처는 {', '.join(SOURCES)} 중 하나입니다.")
    if out["min_precision"] not in PRECISION:
        raise ValueError(f"날짜 정밀도는 {', '.join(PRECISION)} 중 하나입니다.")
    if out["reroot"] not in REROOTS:
        raise ValueError(f"루트는 {', '.join(REROOTS)} 중 하나입니다.")
    try:
        out["clock_filter"] = float(s["clock_filter"])
        out["min_group"] = int(s["min_group"])
    except (TypeError, ValueError):
        raise ValueError("이상치 기준과 최소 그룹 크기는 숫자여야 합니다.")
    if not 0 <= out["clock_filter"] <= 20 or out["min_group"] < 2:
        raise ValueError("이상치 기준은 0~20(0이면 끔), 최소 그룹 크기는 2 이상이어야 합니다.")
    if out["coalescent"] not in ("", "opt", "skyline") and not re.fullmatch(r"\d+(\.\d+)?", out["coalescent"]):
        raise ValueError("coalescent는 비움, opt, skyline 또는 숫자(Tc)입니다.")
    return out


def settings_sha(s: dict) -> str:
    return hashlib.sha256(json.dumps(s, sort_keys=True).encode()).hexdigest()


# ------------------------------------------------------------------------------------------------ inputs: the tree and the dates
def decimal_date(d: date) -> float:
    y0 = date(d.year, 1, 1)
    return d.year + (d - y0).days / ((date(d.year + 1, 1, 1) - y0).days)


def _treetime_date(text: str, min_precision: str) -> tuple[str, str] | None:
    """(date text for TreeTime, precision) of a curated collection date, or None when it is not usable at this precision. A month is the range of that month, a year the range of that year."""
    kind, _week, month, shown = lineages.parse_date(text)
    if kind == "day":
        return shown, "day"
    if kind == "month" and min_precision in ("month", "year"):
        y, m = int(month[:4]), int(month[5:])
        last = (date(y + (m == 12), 1 if m == 12 else m + 1, 1) - date(y, m, 1)).days
        return f"[{y}-{m:02d}-01:{y}-{m:02d}-{last:02d}]", "month"
    if min_precision == "year":
        mm = re.match(r"^(\d{4})$", (text or "").strip())
        if mm:
            return f"[{mm[1]}-01-01:{mm[1]}-12-31]", "year"
    return None


def given_decimal(text: str) -> float:
    """The decimal year of a date as given to TreeTime: a day, or the middle of a range `[a:b]` (never what the regression infers for a tip it left out)."""
    if text.startswith("["):
        a, b = text[1:-1].split(":")
        da, db = date.fromisoformat(a), date.fromisoformat(b)
        return (decimal_date(da) + decimal_date(db)) / 2
    return decimal_date(date.fromisoformat(text))


def _squirrel_tree(src: dict) -> tuple[str, int]:
    d = Path(src["run_dir"]) / "squirrel"
    f = d / "alignment.tree"
    text = f.read_text(encoding="utf-8")
    m = re.search(r"tree\s+\S+\s*=\s*(?:\[&R\]\s*)?(\(.*?;)\s*(?:end;|$)", text, re.S)
    if not m:
        raise ValueError("Squirrel 계통수 파일을 읽지 못했습니다.")
    aln = d / "alignment.aln.fasta.gz"
    length = len(next(aligned.read_alignment(aln))[1]) if aln.exists() else 0
    return m.group(1), length


def _nextstrain_tree(src: dict) -> tuple[str, int]:
    f = nextstrain.runs_dir().parent / "auspice" / src["auspice_file"]
    tree = json.loads(f.read_text(encoding="utf-8"))["tree"]
    ref = nextstrain.workflow_dir() / "defaults" / (src.get("build") or "") / "reference.fasta"
    if not ref.exists():
        ref = nextstrain.workflow_dir() / "defaults" / "reference.fasta"
    from . import seqio
    length = len(seqio.read_fasta(ref)[1]) if ref.exists() else 0
    if not length:
        raise ValueError("Nextstrain 트리의 길이 단위를 정할 참조 서열을 찾지 못했습니다.")

    def nw(n, parent_div):
        div = (n.get("node_attrs") or {}).get("div", parent_div) or 0
        bl = max(0.0, (div - parent_div) / length)
        if "children" in n:
            return "(" + ",".join(nw(c, div) for c in n["children"]) + f"):{bl:.10f}"
        return f"{n['name']}:{bl:.10f}"
    return nw(tree, 0) + ";", length


def find_source(ds: dict, source: str) -> dict | None:
    for r in nextstrain.list_runs(ds["id"]):
        if r["status"] != "done":
            continue
        if source == "squirrel" and r["tool"] == "squirrel":
            try:
                if json.loads(r.get("params_yaml") or "{}").get("run_phylo") is False:
                    continue
            except ValueError:
                pass
            if (Path(r["run_dir"]) / "squirrel" / "alignment.tree").exists():
                return r
        if source == "nextstrain" and r["tool"] == "nextstrain" and r.get("auspice_file"):
            return r
    return None


def tip_names(newick: str) -> list[str]:
    return re.findall(r"[(,]\s*([^(),:;\s]+)\s*:", newick)


# ------------------------------------------------------------------------------------------------ running TreeTime
def _tt(args: list[str], cwd: Path, log_path: Path, timeout: int = 7200) -> None:
    pre = env_prefix()
    cmd = [str(pre / "bin" / "treetime")] + args
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write("$ " + " ".join(cmd) + "\n")
        r = subprocess.run(cmd, cwd=str(cwd), env=nextstrain.tool_env(pre), stdout=lf, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    if r.returncode != 0:
        tail = "".join(open(log_path, encoding="utf-8").readlines()[-12:])
        raise RuntimeError(f"TreeTime 실패 (종료 코드 {r.returncode}): {tail.strip()[-500:]}")


def _num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def parse_clock_txt(path: Path) -> dict:
    t = path.read_text(encoding="utf-8") if path.exists() else ""
    rate = re.search(r"--rate:\s*([0-9.eE+-]+)(?:\s*\+/-\s*([0-9.eE+-]+))?", t)
    r2 = re.search(r"--r\^2:\s*([0-9.eE+-]+)", t)
    return {"rate": _num(rate[1]) if rate else None, "rate_sd": _num(rate[2]) if rate and rate[2] else None, "r2": _num(r2[1]) if r2 else None}


def parse_rtt(path: Path, outliers: Path) -> tuple[list[dict], set[str]]:
    """(tips [{name, date, rtt, dev}], outlier names) from `treetime clock`'s rtt.csv and outliers.tsv; internal nodes are left out."""
    pts = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#") or line.startswith("name"):
            continue
        p = [x.strip() for x in line.split(",")]
        if len(p) >= 4 and not p[0].startswith("NODE_"):
            pts.append({"name": p[0], "date": _num(p[1]), "rtt": _num(p[2]), "dev": _num(p[3])})
    out = set()
    if outliers.exists():
        for line in outliers.read_text(encoding="utf-8").splitlines()[1:]:
            if line.strip():
                out.add(line.split("\t")[0])
    return pts, out


def _walk(node: dict, parent=None, out=None):
    out = out if out is not None else []
    out.append((node, parent))
    for c in node.get("children", []):
        _walk(c, node, out)
    return out


def group_tmrca(tree: dict, groups: dict[str, set[str]], min_size: int) -> list[dict]:
    """The MRCA date (with its interval) of each group of genomes in a TreeTime auspice tree. A group that is not monophyletic gets the MRCA of all its members, which then also holds
    other genomes: `extra` says how many. Groups with fewer than `min_size` genomes in the tree are skipped."""
    nodes = _walk(tree)
    tips_of: dict[int, set[str]] = {}

    def collect(n):
        s = {n["name"]} if "children" not in n else set().union(*(collect(c) for c in n["children"]))
        tips_of[id(n)] = s
        return s
    collect(tree)
    out = []
    for g, members in groups.items():
        m = members & tips_of[id(tree)]
        if len(m) < min_size:
            continue
        best = min((n for n, _p in nodes if m <= tips_of[id(n)]), key=lambda n: len(tips_of[id(n)]))
        nd = (best.get("node_attrs") or {}).get("num_date") or {}
        conf = nd.get("confidence") or [None, None]
        out.append({"group": g, "n": len(m), "extra": len(tips_of[id(best)]) - len(m), "tmrca": nd.get("value"), "lo": conf[0], "hi": conf[1]})
    return sorted(out, key=lambda x: (x["tmrca"] is None, x["tmrca"] or 0))


def warnings(s: dict) -> list[str]:
    """What a reader must keep in mind before using the numbers."""
    w = []
    r2 = (s.get("clock") or {}).get("r2")
    if r2 is not None and r2 < 0.3:
        w.append(f"시계 신호가 약합니다(R² {r2:.2f}). 점추정의 날짜가 있어도 신뢰구간이 매우 넓고, 속도와 루트 날짜는 표본 수집 기간 밖으로 외삽한 값입니다.")
    span = s.get("span_years")
    if span is not None and span < 5:
        w.append(f"표본 수집 기간이 {span:.1f}년뿐입니다. 이 기간에 쌓인 변이가 적으면 속도와 tMRCA를 정밀하게 정할 수 없습니다.")
    root = (s.get("time") or {}).get("root") or {}
    if root.get("lo") is not None and root.get("hi") is not None and root["hi"] - root["lo"] > 2:
        w.append(f"루트 날짜의 90% 구간이 {root['hi'] - root['lo']:.1f}년으로 넓습니다. 이 구간 전체가 가능한 값입니다.")
    if s.get("n_outliers", 0) > 0.03 * max(1, s.get("n_dated", 1)):
        w.append(f"이상치 {s['n_outliers']}개({100 * s['n_outliers'] / max(1, s['n_dated']):.0f}%)가 시계에서 크게 벗어납니다 — 재조합·오염·날짜 오류·APOBEC3 폭발 같은 비시계적 변화일 수 있습니다(제외하고 회귀했습니다).")
    p = s.get("precision") or {}
    if p.get("month") and p["month"] > 0.2 * max(1, s.get("n_dated", 1)):
        w.append(f"날짜가 월 단위인 유전체가 {p['month']}개({100 * p['month'] / max(1, s['n_dated']):.0f}%)입니다. 범위로 넣어 불확실성에 반영했습니다.")
    if s.get("n_dropped_no_date"):
        w.append(f"날짜가 없거나 쓸 수 없는 유전체 {s['n_dropped_no_date']}개는 분석에서 뺐습니다.")
    w.append("이 값은 트리의 루트가 곧 공통 조상이라고 가정합니다. 루트가 잘못 잡히면(외군이 없으면 흔함) 속도와 날짜가 모두 틀어집니다.")
    return w


# ------------------------------------------------------------------------------------------------ the run
def run(ds: dict, settings: dict | None = None, force: bool = False, note: str = "", progress=None) -> dict:
    """Stage 1 + stage 2 on the dataset's tree. An identical finished run (same tree, dates, settings, TreeTime) is returned unless force."""
    p = progress or (lambda *a: None)
    st = status()
    if not st["ok"]:
        raise RuntimeError(f"TreeTime이 없습니다: {st['how']}")
    s = parse_settings(settings)
    src = find_source(ds, s["source"])
    if not src:
        raise ValueError("이 데이터셋에 " + ("계통수까지 만든 Squirrel 분석" if s["source"] == "squirrel" else "Nextstrain 분석") + " 결과가 없습니다. 먼저 그 분석을 실행하세요.")
    members = nextstrain.run_members(src)
    msha = datasets.members_hash(members)
    csha = nextstrain.content_sha(ds)
    sha = settings_sha({**s, "source_run": src["id"]})
    if not force:
        r = runs.reusable(ds, "clock", msha, csha, sha, st["version"])
        if r:
            p(100, f"같은 트리·날짜·설정으로 이미 분석했습니다 ({r['name']})")
            return {**r, "reused": True}
    name, rdir, now = runs.new_run_dir(ds, "clock", "clock")
    logp = rdir / "run.log"
    rid = runs.start(ds, "clock", name, now, ds.get("build") or "", {**s, "source_run": src["id"]}, sha, members, msha, csha, st["version"], rdir, logp, note, {"treetime": st["version"], "source_run": src["id"]})
    try:
        p(3, "트리와 날짜를 준비하는 중…")
        newick, seqlen = (_squirrel_tree if s["source"] == "squirrel" else _nextstrain_tree)(src)
        names = tip_names(newick)
        gmeta = {g["id"]: g for g in lineages._genomes(src)}
        dated, prec, dropped = [], Counter(), 0
        for n in names:
            g = gmeta.get(n)
            d = _treetime_date(g["date"], s["min_precision"]) if g else None
            if d:
                dated.append((n, d[0]))
                prec[d[1]] += 1
            else:
                dropped += 1
        if len(dated) < 10:
            raise ValueError(f"날짜가 있는 유전체가 {len(dated)}개뿐이라 분자시계를 추정할 수 없습니다(최소 10개).")
        (rdir / "tree.nwk").write_text(newick, encoding="utf-8")
        with open(rdir / "dates.csv", "w", encoding="utf-8") as f:
            f.write("name,date\n" + "".join(f"{n},{d}\n" for n, d in dated))
        common = ["--tree", "tree.nwk", "--dates", "dates.csv", "--sequence-length", str(seqlen), "--clock-filter", str(s["clock_filter"])]
        common += ["--keep-root"] if s["reroot"] == "keep" else ["--reroot", s["reroot"]]
        p(8, "1단계: 루트-투-팁 회귀")
        _tt(["clock"] + common + ["--outdir", "rtt"], rdir, logp)                      # stage 1 is TempEst's plain least-squares regression (tips taken as independent)
        clock = parse_clock_txt(rdir / "rtt" / "molecular_clock.txt")
        pts, outl = parse_rtt(rdir / "rtt" / "rtt.csv", rdir / "rtt" / "outliers.tsv")
        root_date_rtt = None
        mm = re.search(r"root-date:\s*([0-9.]+)", logp.read_text(encoding="utf-8"))
        if mm:
            root_date_rtt = float(mm[1])
        p(25, "2단계: 시간 보정 트리와 신뢰구간 (몇 분)")
        t_args = ["--branch-length-mode", "input", "--confidence", "--time-marginal", "assign"] + (["--covariation"] if s["covariation"] else []) + (["--coalescent", s["coalescent"]] if s["coalescent"] else [])
        _tt(common + t_args + ["--outdir", "time"], rdir, logp)             # stage 2 can take the shared ancestry of the tips into account (covariation)
        tc = parse_clock_txt(rdir / "time" / "molecular_clock.txt")
        tr = json.loads((rdir / "time" / "auspice_tree.json").read_text(encoding="utf-8"))["tree"]
        nd = (tr.get("node_attrs") or {}).get("num_date") or {}
        conf = nd.get("confidence") or [None, None]
        groups: dict[str, set] = {}
        for n, _d in dated:
            g = gmeta[n]
            for kind, val in (("clade", g["clade"]), ("lineage", g["lineage"])):
                if val:
                    groups.setdefault(f"{kind} {val}", set()).add(n)
        tmrca = group_tmrca(tr, groups, s["min_group"])
        point = {x["name"]: x for x in pts}
        tips_out = []
        for n, dtxt in dated:
            x = point.get(n)
            if x and x["rtt"] is not None:
                g = gmeta[n]
                tips_out.append({"name": n, "date": given_decimal(dtxt), "rtt": x["rtt"], "dev": x["dev"], "outlier": n in outl, "lineage": g["lineage"], "clade": g["clade"], "country": g["country"],
                                 "precision": "range" if dtxt.startswith("[") else "day"})
        dates = [t["date"] for t in tips_out]
        summary = {"settings": s, "source_run": src["id"], "source_tool": src["tool"], "sequence_length": seqlen, "n_tips": len(names), "n_dated": len(dated), "n_dropped_no_date": dropped,
                   "precision": dict(prec), "n_outliers": len(outl), "clock": {**clock, "root_date": root_date_rtt, "per_genome_per_year": clock["rate"] * seqlen if clock["rate"] else None},
                   "time": {"root": {"date": nd.get("value"), "lo": conf[0], "hi": conf[1]}, "rate": tc["rate"], "rate_sd": tc["rate_sd"], "r2": tc["r2"], "groups": tmrca},
                   "span_years": (max(dates) - min(dates)) if dates else None, "first_sample": min(dates) if dates else None, "last_sample": max(dates) if dates else None, "tips": tips_out,
                   "treetime": st["version"]}
        summary["warnings"] = warnings(summary)
        (rdir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False), encoding="utf-8")
        runs.finish(rid, "clock", name, "done", "", "summary.json")
        p(100, "완료")
    except Exception as e:
        runs.finish(rid, "clock", name, "error", str(e)[:400])
        raise
    return nextstrain.get_run(rid)


def summary(run: dict) -> dict:
    f = Path(run["run_dir"] or "") / "summary.json"
    try:
        return json.loads(f.read_text(encoding="utf-8")) if f.exists() else {}
    except ValueError:
        return {}


def result_files(run: dict) -> list[dict]:
    base = Path(run["run_dir"] or "")
    out = []
    for rel in ("summary.json", "dates.csv", "tree.nwk", "rtt/rtt.csv", "rtt/outliers.tsv", "time/dates.tsv", "time/timetree.nexus", "time/auspice_tree.json", "time/molecular_clock.txt"):
        f = base / rel
        if f.exists():
            out.append({"name": rel, "path": rel, "size": f.stat().st_size})
    return out
