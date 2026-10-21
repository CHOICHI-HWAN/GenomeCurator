"""Nextclade's mutation calls for the genomes of a mutation-profile run, in the run's own reference coordinates, so they can be set beside the run's calls and snippy / lofreq.

Nextclade reports substitutions against *its* dataset reference (for mpox NC_063383.1), the run may use another reference (KJ642613.1 for clade Ib). To compare like with like, a genome's base at
every site where the two references differ is worked out too: a site Nextclade lists as a substitution has the substituted base, a site it reports as missing / N / deleted is unknown (no call, no
absence), any other site has the dataset reference's base. A call is then "the genome's base differs from the run's reference". Only aligned sites (minimap2 liftover of the two references) are
used; the run's ITR window is left out like everywhere else. The result is `mutprofile/nextclade_calls.tsv` (genome, pos, ref, alt) and `mutprofile/nextclade.genomes.txt` (the genomes judged)."""
import csv
import re
import subprocess
import tempfile
from pathlib import Path

from . import config, datasets, lineages, mutprofile, nextstrain, qc, seqio, squirrel
from .db import connect


def _file(run: dict) -> Path:
    return Path(run["run_dir"]) / "mutprofile" / "nextclade_calls.tsv"


def available(run: dict) -> bool:
    return _file(run).exists()


def dataset_of(run: dict) -> Path | None:
    ds = datasets.get(run["dataset_id"])
    with connect() as con:
        w = con.execute("SELECT taxon FROM watchlist WHERE id=?", (ds["watch_id"],)).fetchone()
    path = qc.dataset_for(w["taxon"] if w else None)
    d = qc.dataset_dir(path) if path else None
    return d if d and (d / "reference.fasta").exists() else None


def status(run: dict | None = None) -> dict:
    b = qc.nextclade_bin()
    d = dataset_of(run) if run else None
    return {"ok": bool(b and (d or not run)), "bin": bool(b), "dataset": d.name if d else "", "version": qc.nextclade_version() or ""}


def _ranges(text: str) -> list[tuple[int, int]]:
    """'10-20,300' (or 'N:5-9,R:12' for the non-ACGT column) -> [(10, 20), (300, 300), ...] in dataset-reference coordinates, 1-based inclusive."""
    out = []
    for tok in re.split(r"[,;]", text or ""):
        m = re.search(r"(\d+)(?:-(\d+))?\s*$", tok.strip())
        if m:
            a = int(m[1])
            out.append((a, int(m[2]) if m[2] else a))
    return out


def parse_calls(tsv: Path, nc_seq: str, run_seq: str, to_nc: dict[int, int], window=None) -> dict[str, set]:
    """{genome: {(pos in the run reference, ref, alt)}} from a Nextclade TSV (columns seqName, substitutions, deletions, missing, nonACGTs)."""
    nc_to_run = {q: p for p, q in to_nc.items()}
    refdiff = [p for p, q in to_nc.items() if 0 < p <= len(run_seq) and 0 < q <= len(nc_seq) and run_seq[p - 1].upper() != nc_seq[q - 1].upper() and run_seq[p - 1].upper() in "ACGT" and nc_seq[q - 1].upper() in "ACGT"]
    out: dict[str, set] = {}
    with open(tsv, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["seqName"].split()[0]
            if (r.get("errors") or "").strip():
                continue
            unknown = _ranges(r.get("missing")) + _ranges(r.get("deletions")) + _ranges(r.get("nonACGTs"))
            subs = {}
            for tok in (r.get("substitutions") or "").split(","):
                m = re.fullmatch(r"([ACGT])(\d+)([ACGT])", tok.strip())
                if m:
                    subs[int(m[2])] = m[3]
            calls = set()
            for q, base in subs.items():
                p = nc_to_run.get(q)
                if p and 0 < p <= len(run_seq) and run_seq[p - 1].upper() in "ACGT" and run_seq[p - 1].upper() != base and (window is None or window(p)):
                    calls.add((p, run_seq[p - 1].upper(), base))
            for p in refdiff:                                     # the genome carries the dataset reference's base where Nextclade lists nothing
                q = to_nc[p]
                if q in subs or any(a <= q <= b for a, b in unknown) or (window is not None and not window(p)):
                    continue
                calls.add((p, run_seq[p - 1].upper(), nc_seq[q - 1].upper()))
            out[g] = calls
    return out


def run_for_mutprofile(run: dict, progress=None) -> dict:
    b = qc.nextclade_bin()
    d = dataset_of(run)
    if not b or not d:
        raise RuntimeError("Nextclade 또는 이 바이러스의 Nextclade 데이터셋이 없습니다 (QC 화면에서 데이터셋을 받으세요).")
    ref = mutprofile.summary(run).get("reference", "")
    fa, _gff = mutprofile.reference_by_name(ref)
    run_seq = mutprofile.read_fasta(fa)[1]
    nc_seq = mutprofile.read_fasta(d / "reference.fasta")[1]
    to_nc = seqio.liftover(fa, d / "reference.fasta", mutprofile._bin("minimap2"), mutprofile.tool_env()) if Path(fa).resolve() != (d / "reference.fasta").resolve() else {i: i for i in range(1, len(run_seq) + 1)}
    ds = datasets.get(run["dataset_id"])
    members = set(nextstrain.run_members(run))
    win = mutprofile.run_itr(run)
    inwin = (lambda p: win[0] <= p <= win[1]) if win else None
    with tempfile.TemporaryDirectory(prefix="gcncall-") as td:
        inp, tsv = Path(td) / "in.fasta", Path(td) / "out.tsv"
        squirrel.write_input(ds, inp)
        keep = Path(td) / "keep.fasta"
        n = 0
        with open(keep, "w", encoding="utf-8") as w:
            for name, s in qc._read_fasta(inp):
                acc = name.split()[0]
                if acc in members:
                    w.write(f">{acc}\n{s}\n")
                    n += 1
        if progress:
            progress(5, f"Nextclade {n}건 실행")
        r = subprocess.run([b, "run", "--input-dataset", str(d), "--output-tsv", str(tsv), str(keep)], capture_output=True, text=True, timeout=7200)
        if r.returncode != 0 or not tsv.exists():
            raise RuntimeError(f"Nextclade 실행 실패: {(r.stderr or r.stdout).strip()[-400:]}")
        calls = parse_calls(tsv, nc_seq, run_seq, to_nc, inwin)
    out = _file(run)
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("genome\tpos\tref\talt\n")
        for g in sorted(calls):
            for v in sorted(calls[g]):
                f.write(f"{g}\t{v[0]}\t{v[1]}\t{v[2]}\n")
    out.with_name("nextclade.genomes.txt").write_text("\n".join(sorted(calls)) + "\n", encoding="utf-8")
    return {"genomes": len(calls), "calls": sum(len(v) for v in calls.values()), "reference": ref, "dataset": d.name}


def read_calls(run: dict, judged: set) -> dict[str, set]:
    """{genome: {(pos, ref, alt)}} for the genomes of `judged` that Nextclade was run on."""
    f = _file(run)
    have = set(f.with_name("nextclade.genomes.txt").read_text(encoding="utf-8").split())
    out = {g: set() for g in judged & have}
    for r in csv.DictReader(open(f, encoding="utf-8"), delimiter="\t"):
        if r["genome"] in out:
            out[r["genome"]].add((int(r["pos"]), r["ref"], r["alt"]))
    return out
