"""Lineage definitions (github.com/mpxv-lineages/lineage-designation) and how lineages / defining mutations spread over time and countries.

The definitions repository has one yaml per lineage (`definitions/<clade>/<outbreak>/<lineage>.yml`: parent, designation date, `defining_snps` as position + nucleotide
in **NC_063383.1** coordinates). It is fetched with git into `data/lineage-designation/` on request (a button; the commit is shown). A lineage's full set of defining
SNPs is the chain of its own and its ancestors' (the definitions are stepwise).

Two views of the same genomes: **by lineage** (what the sample's `lineage` column says) and **by defining mutation** (which genomes carry the SNP, whatever lineage they
were given - so the picture stays the same when a lineage is split or renamed later). Carrying a defining SNP is decided on the genome's SNPs against NC_063383.1
(its own pass, kept next to the run), because a run may use another reference (clade I); a base that is N or in a gap counts as *not carrying* (the cautious side).
Dates: the collection date; ISO week when most dates are exact days, otherwise (or when chosen) month - a genome whose date is only a month is never put into a week."""
import csv
import re
import subprocess
from collections import Counter, defaultdict
from datetime import date, timedelta
from pathlib import Path

from . import config, datasets, export, geo, log, mutprofile, nextstrain, seqio, squirrel
from .db import audit, connect

_log = log.get("lineages")
REPO = "https://github.com/mpxv-lineages/lineage-designation"
NC = "NC_063383.1"


def repo_dir() -> Path:
    import os
    return Path(os.environ.get("GC_LINEAGE_DIR") or (config.DATA_DIR / "lineage-designation")).expanduser()


def commit() -> str:
    d = repo_dir()
    git = __import__("shutil").which("git")
    if not git or not (d / ".git").exists():
        return ""
    try:
        return subprocess.run([git, "-C", str(d), "rev-parse", "--short", "HEAD"], capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        return ""


def status() -> dict:
    d = repo_dir()
    n = len(list((d / "definitions").rglob("*.yml"))) if (d / "definitions").exists() else 0
    return {"ok": n > 0, "dir": str(d), "commit": commit(), "files": n, "repo": REPO}


def clone(progress=None) -> dict:
    """git clone --depth 1 of the definitions repository (or a fast-forward pull)."""
    import shutil
    p = progress or (lambda *a: None)
    d = repo_dir()
    git = shutil.which("git") or str((nextstrain.env_prefix() or Path("/x")) / "bin" / "git")
    if (d / ".git").exists():
        p(20, "정의 저장소를 업데이트하는 중…")
        r = subprocess.run([git, "-C", str(d), "pull", "--ff-only"], capture_output=True, text=True, timeout=600)
    else:
        if d.exists() and any(d.iterdir()):
            raise RuntimeError(f"{d} 이(가) 비어 있지 않습니다. 비우거나 GC_LINEAGE_DIR로 다른 위치를 정하세요.")
        d.parent.mkdir(parents=True, exist_ok=True)
        p(20, f"{REPO} 를 받는 중…")
        r = subprocess.run([git, "clone", "--depth", "1", REPO, str(d)], capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        raise RuntimeError(f"git 실패: {(r.stderr or r.stdout).strip()[-300:]}")
    _DEF.clear()
    con = connect()
    try:
        audit(con, "lineage_definitions", f"{REPO} @ {commit()} -> {d}")
        con.commit()
    finally:
        con.close()
    p(100, f"정의 {commit()}")
    return status()


# ------------------------------------------------------------------------------------------------ definitions
_DEF: dict = {}


def definitions() -> dict:
    """{outbreak: {lineage: {name, parent, clade, designation_date, snps: [(pos, nucleotide)]}}}. The outbreak's own yaml (`<clade>/<outbreak>.yml`) is stored as lineage ''."""
    key = str(repo_dir())
    if key in _DEF:
        return _DEF[key]
    import yaml
    out: dict = defaultdict(dict)
    base = repo_dir() / "definitions"
    for f in sorted(base.rglob("*.yml")) if base.exists() else []:
        try:
            y = yaml.safe_load(f.read_text(encoding="utf-8")) or {}
        except Exception:
            continue
        rel = f.relative_to(base).parts
        if len(rel) == 3:                                    # <clade>/<outbreak>/<lineage>.yml
            clade, outbreak, name = rel[0], rel[1], y.get("name") or f.stem
        elif len(rel) == 2:                                  # <clade>/<outbreak>.yml
            clade, outbreak, name = rel[0], f.stem, ""
        else:
            continue
        snps = [(int(s["pos"]), str(s["nucleotide"]).upper()) for s in (y.get("defining_snps") or []) if isinstance(s, dict) and "pos" in s and "nucleotide" in s]
        out[outbreak][name] = {"name": name, "parent": y.get("parent") or "", "clade": clade, "designation_date": str(y.get("designation_date") or ""), "snps": snps}
    _DEF[key] = dict(out)
    return _DEF[key]


def chain(outbreak: str, lineage: str) -> list[dict]:
    """The lineage and its ancestors, root first (the outbreak's own definition first when it has one)."""
    d, seen, out, cur = definitions().get(outbreak, {}), set(), [], lineage
    while cur and cur in d and cur not in seen:
        seen.add(cur)
        out.append(d[cur])
        cur = d[cur]["parent"]
    if "" in d:
        out.append(d[""])
    return out[::-1]


def defining_snps(outbreak: str, lineage: str) -> list[tuple[int, str, str]]:
    """(position in NC_063383.1, nucleotide, lineage that introduced it) for the whole chain."""
    seen, out = set(), []
    for c in chain(outbreak, lineage):
        for pos, nuc in c["snps"]:
            if (pos, nuc) not in seen:                       # a lineage often repeats its parent's SNPs: keep the first (the one that introduced it)
                seen.add((pos, nuc))
                out.append((pos, nuc, c["name"] or outbreak))
    return out


# ------------------------------------------------------------------------------------------------ coordinates
def nc_reference() -> tuple[str, str, list]:
    base = nextstrain.workflow_dir() / "defaults"
    fa, gff = base / "reference.fasta", base / "genome_annotation.gff3"
    if not fa.exists():
        raise RuntimeError("NC_063383.1 참조를 찾을 수 없습니다. Nextstrain 설정 카드에서 워크플로를 먼저 받아 주세요.")
    name, seq = mutprofile.read_fasta(fa)
    return name, seq, mutprofile.load_cds(gff if gff.exists() else None)


def liftover(ref_fa: Path, nc_fa: Path) -> dict[int, int]:
    """{position in `ref_fa`'s reference: position in NC_063383.1} for the aligned (matching or substituted) columns."""
    if str(ref_fa) == str(nc_fa):
        return None
    return seqio.liftover(ref_fa, nc_fa, mutprofile._bin("minimap2"), mutprofile.tool_env())


def gap_intervals_from_paf(paf_lines, length: int) -> dict[str, list[tuple[int, int]]]:
    """{query: [(first, last)]} positions of the reference (1-based) where the genome has no usable base: N / ambiguity characters, deletions and everything the alignment
    does not cover. A defining SNP at such a position is *unknown*, not absent."""
    cover: dict[str, list[tuple[int, int]]] = defaultdict(list)
    gaps: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for line in paf_lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if not cs:
            continue
        q, t = f[0], int(f[7])
        cover[q].append((t + 1, int(f[8])))
        for tok in mutprofile._CS.findall(cs):
            if tok[0] == ":":
                t += int(tok[1:])
            elif tok[0] == "*":
                if tok[2] not in "acgt":
                    gaps[q].append((t + 1, t + 1))
                t += 1
            elif tok[0] == "-":
                gaps[q].append((t + 1, t + len(tok) - 1))
                t += len(tok) - 1
            elif tok[0] == "~":
                n = int(re.search(r"\d+", tok).group())
                gaps[q].append((t + 1, t + n))
                t += n
    for q, cv in cover.items():
        pos = 1
        for a, b in sorted(cv):
            if a > pos:
                gaps[q].append((pos, a - 1))
            pos = max(pos, b + 1)
        if pos <= length:
            gaps[q].append((pos, length))
    return gaps


def nc_unknown(run: dict) -> dict[str, list[tuple[int, int]]]:
    """{genome id: intervals of NC_063383.1 where the genome has no usable base} (built together with the NC_063383.1 SNPs)."""
    _nc_snps(run)
    f = Path(run["run_dir"]) / "mutprofile" / "nc_gaps.tsv"
    out: dict[str, list] = defaultdict(list)
    if f.exists():
        for r in csv.DictReader(open(f, encoding="utf-8"), delimiter="\t"):
            out[r["ID"]].append((int(r["START"]), int(r["END"])))
    return out


def is_unknown(intervals, pos: int) -> bool:
    return any(a <= pos <= b for a, b in intervals or ())


def _nc_snps(run: dict) -> dict[str, dict[int, str]]:
    """{genome id: {position in NC_063383.1: alt}} for the run's members (own pass against NC_063383.1, kept as snps_nc063383.tsv next to the run, with nc_gaps.tsv: where a genome
    has no usable base)."""
    out = Path(run["run_dir"]) / "mutprofile"
    f = out / "snps_nc063383.tsv"
    sm = mutprofile.summary(run)
    if not (f.exists() and (out / "nc_gaps.tsv").exists()) and sm.get("aligner") == "squirrel" and sm.get("alignment") and Path(sm["alignment"]).exists():
        from . import aligned      # the run's SNVs came from Squirrel's alignment: read this table from the same alignment (no second, different alignment)
        nc_fa = nextstrain.workflow_dir() / "defaults" / "reference.fasta"
        ids = {m.split(":", 1)[1] if m.startswith("private:") else m for m in nextstrain.run_members(run)}
        res = aligned.snvs_vs_target(sm["alignment"], squirrel.reference_fasta(sm["squirrel_reference"]), nc_fa, ids)
        with open(f, "w", encoding="utf-8", newline="") as fo:
            w = csv.writer(fo, delimiter="\t", lineterminator="\n")
            w.writerow(["ID", "POS", "REF", "ALT"])
            w.writerows(res["snvs"])
        with open(out / "nc_gaps.tsv", "w", encoding="utf-8", newline="") as fo:
            w = csv.writer(fo, delimiter="\t", lineterminator="\n")
            w.writerow(["ID", "START", "END"])
            for q_, iv in res["gaps"].items():
                w.writerows((q_, a, b) for a, b in iv)
    if not (f.exists() and (out / "nc_gaps.tsv").exists()):
        ds = datasets.get(run["dataset_id"])
        members = set(nextstrain.run_members(run))
        tmp = Path(run["run_dir"]) / "nc_input.fasta"
        squirrel.write_input(ds, tmp)
        keep = {m for m in members}
        # keep only the members of the run (the dataset may have changed since)
        with open(tmp, encoding="utf-8") as fin, open(str(tmp) + ".sel", "w", encoding="utf-8") as fout:
            on = False
            for line in fin:
                if line.startswith(">"):
                    on = line[1:].split()[0] in keep or ("private:" + line[1:].split()[0]) in keep
                if on:
                    fout.write(line)
        nc_fa = nextstrain.workflow_dir() / "defaults" / "reference.fasta"
        paf = Path(run["run_dir"]) / "nc.paf"
        with open(paf, "w") as pf:
            r = subprocess.run([mutprofile._bin("minimap2"), "-x", "asm10", "-c", "--cs", "--secondary=no", "-t", "4", str(nc_fa), str(tmp) + ".sel"], stdout=pf, stderr=subprocess.PIPE, text=True, env=mutprofile.tool_env())
        if r.returncode != 0:
            raise RuntimeError(f"minimap2 실패: {r.stderr[-300:]}")
        with open(f, "w", encoding="utf-8", newline="") as fo:
            w = csv.writer(fo, delimiter="\t", lineterminator="\n")
            w.writerow(["ID", "POS", "REF", "ALT"])
            w.writerows(mutprofile.snvs_from_paf(open(paf, encoding="utf-8")))
        gaps = gap_intervals_from_paf(open(paf, encoding="utf-8"), len(mutprofile.read_fasta(nc_fa)[1]))
        with open(out / "nc_gaps.tsv", "w", encoding="utf-8", newline="") as fo:
            w = csv.writer(fo, delimiter="\t", lineterminator="\n")
            w.writerow(["ID", "START", "END"])
            for q_, iv in gaps.items():
                w.writerows((q_, a, b) for a, b in sorted(iv))
        for x in (tmp, Path(str(tmp) + ".sel"), paf):
            x.unlink(missing_ok=True)
    d: dict[str, dict[int, str]] = defaultdict(dict)
    for r in csv.DictReader(open(f, encoding="utf-8"), delimiter="\t"):
        d[r["ID"]][int(r["POS"])] = r["ALT"]
    return d


# ------------------------------------------------------------------------------------------------ dates
_DAY = re.compile(r"^(\d{4})-(\d\d)-(\d\d)$")
_MON = re.compile(r"^(\d{4})-(\d\d)(?:-XX)?$")


def parse_date(s: str):
    """-> (kind, week_key or None, month_key or None, shown date text): kind is 'day', 'month' or ''. A date with a month but no valid day is month only."""
    s = (s or "").strip()
    m = _DAY.match(s)
    if m:
        try:
            d = date(int(m[1]), int(m[2]), int(m[3]))
            y, w, _ = d.isocalendar()
            return "day", f"{y}-W{w:02d}", f"{d.year}-{d.month:02d}", s
        except ValueError:
            return "", None, None, ""
    m = _MON.match(s)
    if m and 1 <= int(m[2]) <= 12:
        return "month", None, f"{m[1]}-{m[2]}", f"{m[1]}-{m[2]}"
    return "", None, None, ""


def week_monday(key: str) -> date:
    y, w = key.split("-W")
    return date.fromisocalendar(int(y), int(w), 1)


def _weeks(a: str, b: str) -> list[str]:
    d, end, out = week_monday(a), week_monday(b), []
    while d <= end:
        y, w, _ = d.isocalendar()
        out.append(f"{y}-W{w:02d}")
        d += timedelta(days=7)
    return out


def _months(a: str, b: str) -> list[str]:
    y, m = int(a[:4]), int(a[5:])
    out = []
    while f"{y}-{m:02d}" <= b:
        out.append(f"{y}-{m:02d}")
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return out


# ------------------------------------------------------------------------------------------------ the spread
def _genomes(run: dict) -> list[dict]:
    ds = datasets.get(run["dataset_id"])
    _c, rows, _sel, _w = export.table(ds["watch_id"], "nextstrain", ds)
    allr = list(rows)
    if ds.get("include_private"):
        _pc, prow, _f = export.private_table(ds["watch_id"], ds)
        allr += prow
    mem = set(nextstrain.run_members(run))
    return [{"id": r["PPX_accession"], "lineage": r.get("lineage") or "", "outbreak": r.get("outbreak") or "", "clade": r.get("clade") or "", "country": r.get("country") or "",
             "date": r.get("date") or ""} for r in allr if not mem or r["PPX_accession"] in mem or ("private:" + r["PPX_accession"]) in mem]


def spread(run: dict, view: str = "lineage", bin_: str = "auto", outbreak: str = "", min_n: int = 2, lineages=(), countries=()) -> dict:
    """The data of the Gantt-like picture. Rows are lineages (view 'lineage') or defining SNPs (view 'mutation'); each row has per-country sub-rows with counts per
    period, the first / last collection date and the number of genomes. `bin_`: 'week', 'month' or 'auto' (week when the exact-day genomes give at least two per occupied
    country-week on average)."""
    gs = _genomes(run)
    ob_counts = Counter(g["outbreak"] for g in gs if g["outbreak"])
    outbreak = outbreak or (ob_counts.most_common(1)[0][0] if ob_counts else "")
    gs = [g for g in gs if g["outbreak"] == outbreak]
    if lineages:
        gs = [g for g in gs if g["lineage"] in set(lineages)]
    if countries:
        gs = [g for g in gs if g["country"] in set(countries)]
    defs = definitions().get(outbreak, {})
    for g in gs:
        g["kind"], g["week"], g["month"], g["shown"] = parse_date(g["date"])
    dated = [g for g in gs if g["kind"]]
    rows_def: list[tuple[str, dict]] = []                     # (key, info) in the order shown
    info: dict[str, dict] = {}
    if view == "mutation":
        if not defs:
            raise ValueError("정의 파일이 없습니다. 먼저 '정의 받기'를 눌러 lineage-designation을 받아 주세요.")
        ncname, ncseq, cds = nc_reference()
        snps = _nc_snps(run)
        used = sorted({g["lineage"] for g in gs if g["lineage"]})
        want: dict[tuple[int, str], str] = {}
        for lin in used:
            for pos, nuc, by in defining_snps(outbreak, lin):
                want.setdefault((pos, nuc), by)
        members_of: dict[tuple[int, str], list] = defaultdict(list)
        for g in dated:
            sn = snps.get(g["id"])
            if not sn:
                continue
            for (pos, nuc) in want:
                if sn.get(pos) == nuc:
                    members_of[(pos, nuc)].append(g)
        for (pos, nuc), by in sorted(want.items(), key=lambda kv: (defs.get(kv[1], {}).get("designation_date", ""), kv[1], kv[0])):
            ref_b = ncseq[pos - 1] if 0 < pos <= len(ncseq) else "?"
            an = mutprofile.annotate(pos, ref_b, nuc, cds, ncseq) if ref_b in "ACGT" else {"gene": "", "aa_change": "", "effect": ""}
            key = f"{pos}{nuc}"
            info[key] = {"key": key, "label": f"{ref_b}{pos}{nuc}", "by": by or outbreak, "gene": an["gene"], "aa": an["aa_change"], "effect": an["effect"], "snps": []}
            rows_def.append((key, info[key]))
        groups = {k: members_of[(int(re.match(r"\d+", k)[0]), k[-1])] for k, _ in rows_def}
    else:
        by_lin: dict[str, list] = defaultdict(list)
        for g in dated:
            if g["lineage"]:
                by_lin[g["lineage"]].append(g)
        for lin in by_lin:
            d = defs.get(lin, {})
            sn = [(pos, nuc) for pos, nuc, _by in defining_snps(outbreak, lin)] if d else []
            own = [f"{pos}{nuc}" for pos, nuc in d.get("snps", [])]
            info[lin] = {"key": lin, "label": lin, "by": d.get("parent", ""), "designated": d.get("designation_date", ""), "own": own, "snps": [f"{p}{n}" for p, n in sn]}
        rows_def = [(k, info[k]) for k in sorted(by_lin, key=lambda k: min((g["shown"] for g in by_lin[k] if g["kind"] == "day"), default="9"))]
        groups = by_lin
    rows_def = [(k, i) for k, i in rows_def if len(groups.get(k, [])) >= min_n]
    # choose the period
    if bin_ not in ("week", "month"):
        occ = {(k, g["country"], g["week"]) for k, _ in rows_def for g in groups[k] if g["kind"] == "day"}
        nday = sum(1 for k, _ in rows_def for g in groups[k] if g["kind"] == "day")
        bin_ = "week" if occ and nday / len(occ) >= 2 else "month"
    allg = [g for k, _ in rows_def for g in groups[k]]
    wk = sorted({g["week"] for g in allg if g["week"]})
    mo = sorted({g["month"] for g in allg if g["month"]})
    periods = (_weeks(wk[0], wk[-1]) if wk else []) if bin_ == "week" else (_months(mo[0], mo[-1]) if mo else [])
    out_rows = []
    for key, inf in rows_def:
        cs: dict[str, dict] = {}
        for g in groups[key]:
            c = cs.setdefault(g["country"] or "(없음)", {"country": g["country"] or "(없음)", "n": 0, "cells": Counter(), "monthly": Counter(), "dates": []})
            c["n"] += 1
            c["dates"].append(g["shown"])
            if bin_ == "week" and g["kind"] == "day":
                c["cells"][g["week"]] += 1
            elif bin_ == "week":
                c["monthly"][g["month"]] += 1
            else:
                c["cells"][g["month"]] += 1
        crow = []
        for c in sorted(cs.values(), key=lambda c: (min(c["dates"]), c["country"])):
            ds_ = sorted(c["dates"], key=lambda s: (s[:7], s))
            pt = geo.resolve(c["country"])
            crow.append({"country": c["country"], "iso": pt["iso"] if pt else "", "lon": pt["lon"] if pt else None, "lat": pt["lat"] if pt else None, "n": c["n"], "first": ds_[0], "last": ds_[-1], "cells": dict(c["cells"]), "monthly": dict(c["monthly"])})
        alld = sorted((g["shown"] for g in groups[key]), key=lambda s: (s[:7], s))
        out_rows.append({**inf, "n": len(groups[key]), "first": alld[0], "last": alld[-1], "countries": crow})
    starts = [week_monday(k).isoformat() for k in periods] if bin_ == "week" else [k + "-01" for k in periods]
    return {"view": view, "bin": bin_, "periods": periods, "starts": starts, "rows": out_rows, "outbreak": outbreak, "outbreaks": sorted(ob_counts), "genomes": len(gs), "dated": len(dated),
            "day_precision": sum(1 for g in dated if g["kind"] == "day"), "month_only": sum(1 for g in dated if g["kind"] == "month"), "no_date": len(gs) - len(dated),
            "definitions": status()["commit"], "reference": mutprofile.summary(run).get("reference", "")}


def defining_in_run_coordinates(run: dict, ref_fa: Path) -> dict[int, int]:
    """{position in NC_063383.1: position in the run's reference}, to draw defining SNPs on the picture of a run whose reference is not NC_063383.1."""
    nc_fa = nextstrain.workflow_dir() / "defaults" / "reference.fasta"
    if Path(ref_fa).resolve() == nc_fa.resolve():
        return None
    return {nc: q for q, nc in liftover(Path(ref_fa), nc_fa).items()}
