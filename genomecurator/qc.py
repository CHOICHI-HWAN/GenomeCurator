"""Sequence-level QC of the genomes in the database: plain sequence statistics for everything, and Nextclade (clade / lineage and its own
quality score) for the viruses that have a Nextclade dataset.

* Results live in `qc_result` (one row per accession and tool: status, score, metrics, tool/dataset versions, the sequence md5 they were
  computed for). They are derived data: delete them and a new run recomputes them.
* QC never edits a curated value. A Nextclade clade/lineage enters curation as one more *source* (`per_source["nextclade"]`, see
  curate.check): it fills a blank, or - when it disagrees with NCBI / BV-BRC / Pathoplexus - shows up in the review screen to be chosen.
* Warnings only: a bad QC status never blocks anything.
* Nextclade datasets are downloaded on demand from Nextstrain (the intended way); they are not part of this repository."""
import csv
import json
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

from . import config, log
from .db import audit, connect

_log = log.get("qc")
RANK = {"good": 0, "mediocre": 1, "bad": 2, "error": 3}
N_WARN, N_BAD = 1.0, 5.0              # % of N in the sequence
IUPAC_WARN = 0.1                      # % of non-ACGTN ambiguity codes
BATCH = 1500                          # sequences per nextclade call (bounded memory, visible progress)

# NCBI species name (lower case, as stored in watchlist.taxon) -> Nextclade dataset. Only datasets that cover a whole species are listed;
# add or override with GC_NEXTCLADE_DATASETS="Species name=dataset/path;Other species=other/path".
DATASETS = {
    "orthopoxvirus monkeypox": "nextstrain/mpox/all-clades",
    "orthoebolavirus zairense": "nextstrain/orthoebolavirus/ebov",
    "orthoebolavirus sudanense": "nextstrain/orthoebolavirus/sudv",
    "orthoebolavirus bundibugyoense": "nextstrain/orthoebolavirus/bdbv",
    "morbillivirus hominis": "nextstrain/measles/genome/WHO-2012",
    "orthoflavivirus denguei": "nextstrain/dengue/all",
    "orthoflavivirus nilense": "nextstrain/wnv/all-lineages",
    "rubivirus rubellae": "nextstrain/rubella/genome",
    "varicellovirus humanalpha3": "nextstrain/herpes/vzv/NC_001348",
    "metapneumovirus hominis": "nextstrain/hmpv/all-clades/NC_039199",
}


def dataset_for(taxon: str | None) -> str | None:
    t = (taxon or "").strip().lower()
    extra = {}
    for part in (os.environ.get("GC_NEXTCLADE_DATASETS") or "").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            extra[k.strip().lower()] = v.strip()
    return extra.get(t) or DATASETS.get(t)


# ------------------------------------------------------------------------------------------------ tool and dataset management
def nextclade_bin() -> str | None:
    return shutil.which("nextclade")


def nextclade_version() -> str | None:
    b = nextclade_bin()
    if not b:
        return None
    try:
        out = subprocess.run([b, "--version"], capture_output=True, text=True, timeout=20).stdout.strip()
        return out.split()[-1] if out else None
    except Exception:
        return None


def dataset_dir(path: str) -> Path:
    return config.DATA_DIR / "nextclade" / path.replace("/", "__")


def local_tag(path: str) -> str | None:
    pj = dataset_dir(path) / "pathogen.json"
    try:
        return json.loads(pj.read_text(encoding="utf-8")).get("version", {}).get("tag")
    except (OSError, ValueError):
        return None


def latest_tag(path: str) -> str | None:
    """The newest dataset version Nextstrain offers (needs internet); None when it cannot be reached."""
    b = nextclade_bin()
    if not b:
        return None
    try:
        r = subprocess.run([b, "dataset", "list", "--name", path, "--json"], capture_output=True, text=True, timeout=90)
        d = json.loads(r.stdout or "[]")
        return d[0]["version"]["tag"] if d else None
    except Exception:
        _log.warning("could not read the latest version of dataset %s", path, exc_info=True)
        return None


def dataset_update(path: str, progress=None) -> dict:
    """Download the newest version of a dataset (into a temp folder first, so a failed download leaves the old one in place)."""
    p = progress or (lambda *a: None)
    b = nextclade_bin()
    if not b:
        raise RuntimeError("Nextclade가 설치되어 있지 않습니다. `mamba install -c bioconda nextclade` (또는 environment.yml로 환경 다시 만들기)")
    before = local_tag(path)
    dest = dataset_dir(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="ncds-", dir=dest.parent))
    try:
        p(20, f"데이터셋 {path} 내려받는 중…")
        r = subprocess.run([b, "dataset", "get", "--name", path, "--output-dir", str(tmp / "d")], capture_output=True, text=True, timeout=900)
        if r.returncode != 0 or not (tmp / "d" / "pathogen.json").exists():
            raise RuntimeError(f"데이터셋을 받지 못했습니다: {(r.stderr or r.stdout).strip()[-300:]}")
        if dest.exists():
            shutil.rmtree(dest)
        shutil.move(str(tmp / "d"), str(dest))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    after = local_tag(path)
    con = connect()
    try:
        audit(con, "qc_dataset", f"{path}: {before or '(없음)'} -> {after}")
        con.commit()
    finally:
        con.close()
    p(100, f"데이터셋 {path} 버전 {after}")
    return {"dataset": path, "before": before, "after": after}


def dataset_status() -> list[dict]:
    """One row per watched virus: which Nextclade dataset applies, which version is on disk, what is the newest (None = unknown)."""
    con = connect()
    try:
        ws = [dict(r) for r in con.execute("SELECT id, name, taxon FROM watchlist ORDER BY name")]
        n = {r[0]: r[1] for r in con.execute("SELECT watch_id, COUNT(*) FROM genome_sequence s JOIN genome g ON g.genome_id=s.genome_id GROUP BY watch_id")}
    finally:
        con.close()
    out, latest_cache = [], {}
    inst = nextclade_bin() is not None
    for w in ws:
        ds = dataset_for(w["taxon"])
        row = {"watch_id": w["id"], "name": w["name"], "taxon": w["taxon"], "dataset": ds, "n_seq": n.get(w["id"], 0), "tool": inst}
        if ds:
            row["local"] = local_tag(ds)
            row["latest"] = latest_cache.setdefault(ds, latest_tag(ds) if inst else None)
            row["state"] = ("missing" if not row["local"] else "unknown" if not row["latest"] else "outdated" if row["latest"] != row["local"] else "current")
        else:
            row["local"] = row["latest"] = None
            row["state"] = "none"
        out.append(row)
    return out


# ------------------------------------------------------------------------------------------------ sequence statistics
def _read_fasta(path: Path):
    """Yield (header, sequence) of a FASTA file (normally one record)."""
    head, buf = None, []
    with open(path, encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith(">"):
                if head is not None:
                    yield head, "".join(buf)
                head, buf = ln[1:].strip(), []
            else:
                buf.append(ln.strip())
    if head is not None:
        yield head, "".join(buf)


def seqstats(path: Path, recorded_length: int | None = None) -> dict:
    """Plain statistics and a verdict (good / mediocre / bad) with the reasons."""
    seq = "".join(s for _, s in _read_fasta(path)).upper()
    n = len(seq)
    if not n:
        return {"status": "bad", "score": None, "metrics": {"length": 0}, "reasons": ["서열이 비어 있습니다"]}
    c = {k: seq.count(k) for k in "ACGTN"}
    other = n - sum(c.values())
    n_pct, iu_pct = 100 * c["N"] / n, 100 * other / n
    gc = 100 * (c["G"] + c["C"]) / max(1, c["A"] + c["C"] + c["G"] + c["T"])
    reasons, level = [], 0
    if n_pct >= N_BAD:
        level = 2
        reasons.append(f"N이 {n_pct:.1f}%")
    elif n_pct >= N_WARN:
        level = max(level, 1)
        reasons.append(f"N이 {n_pct:.1f}%")
    if iu_pct >= IUPAC_WARN:
        level = max(level, 1)
        reasons.append(f"모호 염기 {iu_pct:.2f}%")
    if recorded_length and recorded_length != n:
        level = 2
        reasons.append(f"기록된 길이 {recorded_length:,} ≠ 서열 길이 {n:,}")
    return {"status": ["good", "mediocre", "bad"][level], "score": round(n_pct, 3),
            "metrics": {"length": n, "n": c["N"], "n_pct": round(n_pct, 3), "ambiguous": other, "ambiguous_pct": round(iu_pct, 3), "gc": round(gc, 2)},
            "reasons": reasons}


# ------------------------------------------------------------------------------------------------ nextclade
def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def parse_nextclade_tsv(path: Path) -> dict[str, dict]:
    """{accession: {status, score, clade, lineage, calls, metrics, error}}. The columns between `seqName` and `qc.overallScore`
    are the dataset's own classifications (clade, lineage, outbreak, ...)."""
    out = {}
    with open(path, encoding="utf-8", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        cols = rd.fieldnames or []
        i0, i1 = (cols.index("seqName") + 1 if "seqName" in cols else 0), (cols.index("qc.overallScore") if "qc.overallScore" in cols else 0)
        custom = [c for c in cols[i0:i1] if not c.startswith("mutationPatterns")]
        for r in rd:
            acc = r["seqName"].split()[0]
            err = (r.get("errors") or "").strip()
            calls = {c: r[c] for c in custom if r.get(c)}
            st = (r.get("qc.overallStatus") or "").strip().lower()
            reasons = [f"{c.split('.')[1]} {r[c]}" for c in cols if re.fullmatch(r"qc\.\w+\.status", c) and (r.get(c) or "").lower() in ("mediocre", "bad")]
            out[acc] = {"status": "error" if err else (st if st in RANK else "error"), "score": _num(r.get("qc.overallScore")), "reasons": reasons or ([err] if err else []),
                        "clade": calls.get("clade", ""), "lineage": calls.get("lineage", ""), "calls": calls, "error": err,
                        "metrics": {k: _num(r.get(k)) for k in ("coverage", "totalMissing", "totalSubstitutions", "totalFrameShifts", "alignmentScore")}}
    return out


def _run_nextclade(ds_dir: Path, items: list[tuple[str, Path]]) -> dict[str, dict]:
    b = nextclade_bin()
    with tempfile.TemporaryDirectory(prefix="gcnc-") as td:
        fa, tsv = Path(td) / "in.fa", Path(td) / "out.tsv"
        with open(fa, "w", encoding="utf-8") as w:
            for acc, p in items:
                for k, (_, s) in enumerate(_read_fasta(p)):
                    w.write(f">{acc}\n{s}\n")
                    break          # one record per genome file
        r = subprocess.run([b, "run", "--input-dataset", str(ds_dir), "--output-tsv", str(tsv), str(fa)], capture_output=True, text=True, timeout=3600)
        if r.returncode != 0 or not tsv.exists():
            raise RuntimeError(f"Nextclade 실행 실패: {(r.stderr or r.stdout).strip()[-400:]}")
        return parse_nextclade_tsv(tsv)


# ------------------------------------------------------------------------------------------------ the QC run
def _store(con, acc, tool, watch_id, res, version, dataset, tag, md5):
    con.execute("INSERT OR REPLACE INTO qc_result(accession,tool,watch_id,status,score,clade,lineage,metrics_json,version,dataset,dataset_tag,seq_md5,checked_at)"
                " VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (acc, tool, watch_id, res["status"], res.get("score"), res.get("clade") or None, res.get("lineage") or None,
                 json.dumps({"metrics": res.get("metrics", {}), "reasons": res.get("reasons", []), "calls": res.get("calls", {}), "error": res.get("error", "")},
                            ensure_ascii=False), version, dataset, tag, md5, datetime.now().isoformat(timespec="seconds")))


def run(watch_id: int, progress=None, force: bool = False) -> dict:
    """Sequence stats for every genome of the watch that has a sequence file, and Nextclade when the virus has a dataset.
    Unchanged sequences (same md5, same tool and dataset version) are skipped unless force."""
    p = progress or (lambda *a: None)
    con = connect()
    try:
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        if not w:
            raise ValueError("바이러스를 찾을 수 없습니다.")
        rows = [dict(r) for r in con.execute(
            "SELECT g.accession, g.genome_length, s.md5, s.file_path FROM genome g JOIN genome_sequence s ON s.genome_id=g.genome_id WHERE g.watch_id=?", (watch_id,))]
        have = {(r[0], r[1]): (r[2], r[3], r[4]) for r in con.execute(
            "SELECT accession, tool, seq_md5, version, dataset_tag FROM qc_result WHERE watch_id=?", (watch_id,))}
    finally:
        con.close()
    todo = []
    for r in rows:
        f = config.data_path(r["file_path"]) if r["file_path"] else None
        if f and f.exists():
            todo.append((r, f))
    summary = {"virus": w["name"], "sequences": len(todo), "stats": 0, "nextclade": 0, "skipped": 0, "dataset": None, "notes": [], "no_file": len(rows) - len(todo)}
    if not todo:
        summary["notes"].append("서열 파일이 있는 유전체가 없습니다 (유전체 화면의 Get FASTA로 먼저 받으세요).")
        return summary
    ds = dataset_for(w["taxon"])
    ver = nextclade_version()
    tag = local_tag(ds) if ds else None
    use_nc = bool(ds and ver and tag)
    if ds and not ver:
        summary["notes"].append("이 바이러스는 Nextclade 데이터셋이 있지만 Nextclade가 설치되어 있지 않아 서열 통계만 했습니다.")
    elif ds and not tag:
        summary["notes"].append(f"Nextclade 데이터셋 {ds}이(가) 아직 내려받아지지 않았습니다. QC 탭에서 '데이터셋 업데이트'를 누르세요. 서열 통계만 했습니다.")
    elif not ds:
        summary["notes"].append("이 바이러스에는 연결된 Nextclade 데이터셋이 없어 서열 통계만 했습니다.")
    summary["dataset"] = ds if use_nc else None
    # 1) sequence statistics
    con = connect()
    try:
        for i, (r, f) in enumerate(todo):
            key = (r["accession"], "seqstats")
            if not force and have.get(key, (None,))[0] == r["md5"]:
                summary["skipped"] += 1
                continue
            try:
                res = seqstats(f, r["genome_length"])
            except Exception as e:
                res = {"status": "error", "score": None, "metrics": {}, "reasons": [f"서열을 읽지 못했습니다: {e}"]}
            _store(con, r["accession"], "seqstats", watch_id, res, "builtin", None, None, r["md5"])
            summary["stats"] += 1
            if i % 200 == 0:
                con.commit()
                p(5 + 35 * i / len(todo), f"서열 통계 {i:,}/{len(todo):,}")
        con.commit()
    finally:
        con.close()
    # 2) nextclade
    if use_nc:
        pend = [(r, f) for r, f in todo if force or have.get((r["accession"], "nextclade"), (None, None, None)) != (r["md5"], ver, tag)]
        summary["skipped"] += len(todo) - len(pend)
        for k in range(0, len(pend), BATCH):
            chunk = pend[k:k + BATCH]
            p(40 + 55 * k / max(1, len(pend)), f"Nextclade {k:,}/{len(pend):,} ({ds})")
            res = _run_nextclade(dataset_dir(ds), [(r["accession"], f) for r, f in chunk])
            con = connect()
            try:
                for r, _ in chunk:
                    x = res.get(r["accession"]) or {"status": "error", "score": None, "error": "결과에 없음(정렬 실패)", "metrics": {}}
                    _store(con, r["accession"], "nextclade", watch_id, x, ver, ds, tag, r["md5"])
                    summary["nextclade"] += 1
                con.commit()
            finally:
                con.close()
    con = connect()
    try:
        audit(con, "qc", f"{w['name']}: 통계 {summary['stats']}건, Nextclade {summary['nextclade']}건, 건너뜀 {summary['skipped']}건")
        con.commit()
    finally:
        con.close()
    p(98, "QC 완료")
    _log.info("qc %r: %s", w["name"], {k: v for k, v in summary.items() if k != "notes"})
    return summary


# ------------------------------------------------------------------------------------------------ reading results
def nextclade_source(watch_id: int, con=None) -> dict[str, dict]:
    """{accession: {clade, lineage}} of the current Nextclade results, for use as one more source in curate.check.
    Only results computed for the sequence that is stored now (md5 matches) and not failed."""
    own = con is None
    con = con or connect()
    try:
        out = {}
        for r in con.execute("SELECT q.accession, q.clade, q.lineage FROM qc_result q JOIN genome g ON g.accession=q.accession AND g.watch_id=q.watch_id "
                             "JOIN genome_sequence s ON s.genome_id=g.genome_id AND s.md5=q.seq_md5 "
                             "WHERE q.tool='nextclade' AND q.watch_id=? AND q.status!='error'", (watch_id,)):
            rec = {k: r[k] for k in ("clade", "lineage") if r[k]}
            if rec:
                out[r["accession"]] = rec
        return out
    finally:
        if own:
            con.close()


def summary(watch_id: int | None = None) -> dict:
    """Counts for the QC page: per tool and status, plus the worst genomes."""
    con = connect()
    try:
        q = "SELECT tool, status, COUNT(*) n FROM qc_result" + (" WHERE watch_id=?" if watch_id else "") + " GROUP BY tool, status"
        by = {}
        for r in con.execute(q, (watch_id,) if watch_id else ()):
            by.setdefault(r["tool"], {})[r["status"]] = r["n"]
        w = " AND q.watch_id=?" if watch_id else ""
        worst = [dict(r) for r in con.execute(
            "SELECT q.accession, q.tool, q.status, q.score, q.metrics_json, q.watch_id FROM qc_result q WHERE q.status IN ('bad','error','mediocre')" + w +
            " ORDER BY CASE q.status WHEN 'error' THEN 0 WHEN 'bad' THEN 1 ELSE 2 END, q.score DESC LIMIT 60", (watch_id,) if watch_id else ())]
        for x in worst:
            m = json.loads(x.pop("metrics_json") or "{}")
            x["reasons"] = m.get("reasons") or ([m["error"]] if m.get("error") else [])
        clade = [dict(r) for r in con.execute(
            "SELECT g.accession, g.clade, g.lineage, q.clade nc_clade, q.lineage nc_lineage FROM qc_result q JOIN genome g ON g.accession=q.accession AND g.watch_id=q.watch_id "
            "WHERE q.tool='nextclade' AND q.status!='error'" + w.replace("q.watch_id", "q.watch_id"), (watch_id,) if watch_id else ())]
        diff = [c for c in clade if (c["nc_clade"] and c["clade"] and c["nc_clade"] != c["clade"]) or (c["nc_lineage"] and c["lineage"] and c["nc_lineage"] != c["lineage"])]
        fill = [c for c in clade if (c["nc_clade"] and not c["clade"]) or (c["nc_lineage"] and not c["lineage"])]
        return {"by": by, "worst": worst, "clade_diff": len(diff), "clade_fill": len(fill), "clade_examples": diff[:15], "fill_examples": fill[:10], "total_nextclade": len(clade)}
    finally:
        con.close()


def per_genome(watch_id: int) -> dict[str, dict]:
    """{accession: {qc: worst status over the tools, nc_clade, nc_lineage, nc_status, stats_status}} for the genome tables."""
    con = connect()
    try:
        out: dict[str, dict] = {}
        for r in con.execute("SELECT accession, tool, status, clade, lineage FROM qc_result WHERE watch_id=?", (watch_id,)):
            d = out.setdefault(r["accession"], {"qc": "good"})
            if RANK[r["status"]] > RANK[d["qc"]]:
                d["qc"] = r["status"]
            if r["tool"] == "nextclade":
                d.update(nc_status=r["status"], nc_clade=r["clade"], nc_lineage=r["lineage"])
            else:
                d["stats_status"] = r["status"]
        return out
    finally:
        con.close()
