"""Small sequence helpers with no dependency on the rest of the app, so that the modules that need them (mutation profile, lineages, Squirrel alignments) do not have to import each other."""
import re
import subprocess
from pathlib import Path

CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")      # the pieces of a minimap2 `cs` string


def read_fasta(path: Path) -> tuple[str, str]:
    """(name, upper-case sequence) of the first record of a FASTA file."""
    name, seq = "", []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if name:
                break
            name = line[1:].split()[0]
        else:
            seq.append(line.strip())
    return name, "".join(seq).upper()


def liftover(ref_fa: Path, nc_fa: Path, minimap2: str, env: dict | None = None) -> dict[int, int]:
    """{position in `ref_fa`'s sequence: position in `nc_fa`'s} for the aligned (matching or substituted) columns, by a minimap2 genome-to-genome alignment."""
    cmd = [minimap2, "-x", "asm20", "-c", "--cs", str(nc_fa), str(ref_fa)]
    out = subprocess.run(cmd, capture_output=True, text=True, env=env).stdout
    m: dict[int, int] = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12 or f[4] != "+":
            continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if not cs:
            continue
        q, t = int(f[2]), int(f[7])                          # 0-based starts
        for tok in CS.findall(cs):
            if tok[0] == ":":
                n = int(tok[1:])
                for i in range(n):
                    m[q + i + 1] = t + i + 1
                q += n
                t += n
            elif tok[0] == "*":
                m[q + 1] = t + 1
                q += 1
                t += 1
            elif tok[0] == "+":
                q += len(tok) - 1
            elif tok[0] == "-":
                t += len(tok) - 1
    return m
