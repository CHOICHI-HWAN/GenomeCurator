"""SNVs read from a Squirrel alignment instead of an alignment this app makes itself.

Squirrel (github.com/aineniamh/squirrel) aligns every genome to its own reference for the clade (PX667572 for clade I, NC_063383 for clade II; the output has one row per genome and
one column per reference position, insertions removed, problem sites masked as N). For monkeypox this alignment is the reliable one (Nextclade's has errors for it), so SNVs are called
from it. The columns are put on the reference the analysis uses (a reference-to-reference alignment by minimap2 gives the column map), so positions, genes and lineage definitions keep
meaning what they meant: a genome's base at a column is compared with that reference's base.

A base that is N / ambiguous, and every reference position the alignment does not cover, is *unknown* - never a difference and never "same as the reference"."""
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np

from . import lineages, mutprofile

ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)


def read_alignment(path):
    """(genome id, sequence) of every record; the id is the first word of the header. Plain or gzipped FASTA."""
    op = gzip.open if str(path).endswith(".gz") else open
    name, parts = None, []
    with op(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def _runs(positions: np.ndarray) -> list[tuple[int, int]]:
    """Sorted 1-based positions -> [(first, last)] of consecutive stretches."""
    if len(positions) == 0:
        return []
    cut = np.where(np.diff(positions) != 1)[0]
    starts = np.r_[positions[0], positions[cut + 1]]
    ends = np.r_[positions[cut], positions[-1]]
    return [(int(a), int(b)) for a, b in zip(starts, ends)]


def snvs_vs_target(aln_path, aln_ref_fa, target_fa, ids: set | None = None, progress=None) -> dict:
    """SNVs of every genome of the alignment against the target reference.

    -> {"snvs": [(id, pos, ref, alt)], "gaps": {id: [(first, last)]}, "n": genomes, "covered": share of the target the column map reaches, "length_ok": bool}
    Raises ValueError when a row is not as long as the alignment reference (a different alignment than the one described)."""
    _n1, aseq = mutprofile.read_fasta(Path(aln_ref_fa))
    _n2, tseq = mutprofile.read_fasta(Path(target_fa))
    cmap = lineages.liftover(Path(aln_ref_fa), Path(target_fa)) if str(aln_ref_fa) != str(target_fa) else {i: i for i in range(1, len(tseq) + 1)}
    qs = np.array(sorted(cmap), dtype=np.int64)
    ts = np.array([cmap[q] for q in qs], dtype=np.int64)
    ts_ok = ts <= len(tseq)
    qs, ts = qs[ts_ok], ts[ts_ok]
    tb = np.frombuffer(tseq.upper().encode(), dtype=np.uint8)[ts - 1]
    t_acgt = np.isin(tb, ACGT)
    uncovered = np.setdiff1d(np.arange(1, len(tseq) + 1), ts)
    uncovered_runs = _runs(uncovered)
    snvs: list[tuple[str, int, str, str]] = []
    gaps: dict[str, list] = {}
    n = 0
    for name, seq in read_alignment(aln_path):
        if ids is not None and name not in ids:
            continue
        if len(seq) != len(aseq):
            raise ValueError(f"정렬 행 {name}의 길이({len(seq)})가 Squirrel 참조({len(aseq)})와 다릅니다 — 다른 정렬 파일일 수 있습니다.")
        a = np.frombuffer(seq.upper().encode(), dtype=np.uint8)[qs - 1]
        a_acgt = np.isin(a, ACGT)
        diff = np.where(a_acgt & t_acgt & (a != tb))[0]
        for j in diff:
            snvs.append((name, int(ts[j]), chr(tb[j]), chr(a[j])))
        unknown = ts[~a_acgt]
        gaps[name] = sorted(set(_runs(np.sort(unknown))) | set(uncovered_runs))
        n += 1
        if progress and n % 100 == 0:
            progress(n)
    return {"snvs": snvs, "gaps": gaps, "n": n, "covered": round(len(ts) / max(1, len(tseq)), 4), "length_ok": True}


def reproject(aln_path, aln_ref_fa, target_fa, ids: set | None = None):
    """Yield (id, aligned sequence) of the Squirrel alignment *on the target reference's columns*: the sequence has exactly the target's length, a column holds the genome's base at the
    aligned Squirrel column, and every target position the column map does not reach is N. This is the shape of what Nextclade writes with `--output-fasta` (aligned to the reference,
    insertions relative to the reference removed), so a workflow that expects an alignment to its own reference can take it as it is."""
    _n1, aseq = mutprofile.read_fasta(Path(aln_ref_fa))
    _n2, tseq = mutprofile.read_fasta(Path(target_fa))
    cmap = lineages.liftover(Path(aln_ref_fa), Path(target_fa)) if str(aln_ref_fa) != str(target_fa) else {i: i for i in range(1, len(tseq) + 1)}
    qs = np.array(sorted(cmap), dtype=np.int64)
    ts = np.array([cmap[q] for q in qs], dtype=np.int64)
    keep = ts <= len(tseq)
    qs, ts = qs[keep], ts[keep]
    for name, seq in read_alignment(aln_path):
        if ids is not None and name not in ids:
            continue
        if len(seq) != len(aseq):
            raise ValueError(f"정렬 행 {name}의 길이({len(seq)})가 Squirrel 참조({len(aseq)})와 다릅니다.")
        out = np.full(len(tseq), ord("N"), dtype=np.uint8)
        out[ts - 1] = np.frombuffer(seq.upper().encode(), dtype=np.uint8)[qs - 1]
        yield name, out.tobytes().decode()
