"""A shared cache of the per-genome variant calls (snippy -> lofreq), so that a genome is called once per reference, not once per analysis run.

Calls depend only on the genome's sequence, the reference and the tools' versions - not on the dataset or the settings of a run. The cache therefore keeps, per
`<reference name>-<reference md5>-<snippy version>-<lofreq version>`, one small file per genome sequence (named by the md5 of the upper-cased sequence): the rows
`genome pos ref alt snippy lofreq dp af` as the flow writes them, with the genome id replaced on the way out. A run that meets a genome it has called before copies the rows into its
own `calls/` folder, and Snakemake - finding the result already there - does not call it again. New data therefore costs only its new genomes.

Nothing in the cache is ever trusted blindly: a different reference or tool version is a different cache folder, and a missing or unreadable file simply means "call it again"."""
import hashlib
import re
from pathlib import Path

from . import config


def seq_md5(seq: str) -> str:
    return hashlib.md5(seq.upper().encode()).hexdigest()


def safe_id(gid: str) -> str:
    """The file-name form of a genome id (the same rule as the Snakefile's `SAFE`)."""
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", gid)


def key(ref_name: str, ref_seq: str, snippy: str, lofreq: str) -> str:
    clean = lambda s: re.sub(r"[^A-Za-z0-9_.\-]", "_", s or "unknown")
    return f"{clean(ref_name)}-{hashlib.md5(ref_seq.upper().encode()).hexdigest()[:10]}-snippy{clean(snippy)}-lofreq{clean(lofreq)}"


def folder(k: str) -> Path:
    return config.DATA_DIR / "variantcache" / k


def prefill(records: list[tuple[str, str]], k: str, calls_dir: Path) -> int:
    """Copy the cached calls of the genomes into `calls_dir/<safe id>.tsv` (not over a file that is already there). `records` = [(genome id, sequence)]. -> how many were taken from the cache."""
    calls_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for gid, seq in records:
        dest = calls_dir / f"{safe_id(gid)}.tsv"
        src = folder(k) / f"{seq_md5(seq)}.tsv"
        if dest.exists() or not src.exists():
            continue
        try:
            rows = src.read_text(encoding="utf-8").splitlines()
        except OSError:
            continue
        dest.write_text("".join(f"{gid}\t{r.split(chr(9), 1)[1]}\n" for r in rows if "\t" in r), encoding="utf-8")
        n += 1
    return n


def store(records: list[tuple[str, str]], k: str, calls_dir: Path) -> int:
    """Put the calls of the genomes that a finished flow wrote into the cache (the genome id is dropped from the first column). -> how many new files."""
    d = folder(k)
    n = 0
    for gid, seq in records:
        src = calls_dir / f"{safe_id(gid)}.tsv"
        dest = d / f"{seq_md5(seq)}.tsv"
        if dest.exists() or not src.exists():
            continue
        d.mkdir(parents=True, exist_ok=True)
        text = "".join(f"-\t{r.split(chr(9), 1)[1]}\n" for r in src.read_text(encoding="utf-8").splitlines() if "\t" in r)
        tmp = dest.with_suffix(".tmp")
        tmp.write_text(text, encoding="utf-8")
        tmp.replace(dest)                                              # a half-written file is never seen under its real name
        n += 1
    return n
