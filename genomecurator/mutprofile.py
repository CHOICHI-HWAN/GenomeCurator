"""Mutation profile of a dataset (github.com/insapathogenomics/mutation_profile, GPL-3.0, run as a separate program) with the gene each mutation falls in.

Pipeline (public and private members together, like the other analyses): every genome is aligned to the build's reference (minimap2, `--cs`), its SNVs are read
off the alignment, each SNV is placed in the reference annotation (gene, codon, amino-acid change, synonymous or not), `mutation-profile` adds the flanking
sequence context and the profile (APOBEC3-type GA>AA / TC>TT by default), and the per-mutation and per-gene tables are written. The tools live in the `gc-squirrel`
conda environment (`pip install mutation-profile` is part of squirrel/environment.yml). Runs are kept in `analysis_run` (tool = 'mutprofile') like the others."""
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from urllib.parse import unquote

from . import datasets, export, log, nextstrain, squirrel
from .db import audit, connect

_log = log.get("mutprofile")
BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AMINO[i] for i, (a, b, c) in enumerate((x, y, z) for x in BASES for y in BASES for z in BASES)}
COMP = str.maketrans("ACGT", "TGCA")


def tool_env() -> dict:
    return squirrel.tool_env()


def _bin(name: str) -> str | None:
    pre = squirrel.env_prefix()
    return str(pre / "bin" / name) if pre and (pre / "bin" / name).exists() else None


def status() -> dict:
    ok = bool(_bin("mutation-profile") and _bin("minimap2"))
    ver = ""
    if ok:
        try:
            ver = subprocess.run([_bin("mutation-profile"), "-v"], capture_output=True, text=True, timeout=60, env=tool_env()).stdout.strip().replace("version:", "").strip()
        except Exception:
            pass
    return {"ok": ok, "version": ver, "how": f"mamba env update -f squirrel/environment.yml -n {squirrel.env_name()}   (또는 환경 안에서 pip install --no-deps mutation-profile)"}


# ------------------------------------------------------------------------------------------------ settings
def defaults_for(ds: dict) -> dict:
    return {"before": 5, "after": 5, "profiles": "GA>AA,TC>TT", "reference": "auto"}


def parse_settings(text) -> dict:
    if isinstance(text, dict):
        d = dict(text)
    elif not (text or "").strip():
        return {}
    else:
        try:
            d = json.loads(text)
        except ValueError as e:
            raise ValueError(f"설정(JSON)을 읽을 수 없습니다: {e}") from None
    if not isinstance(d, dict):
        raise ValueError("설정은 JSON 객체여야 합니다.")
    unknown = sorted(set(d) - {"before", "after", "profiles", "reference"})
    if unknown:
        raise ValueError(f"알 수 없는 설정: {', '.join(unknown)}")
    out = {}
    for k in ("before", "after"):
        if str(d.get(k, "")).strip() != "":
            if not str(d[k]).isdigit() or not 0 <= int(d[k]) <= 50:
                raise ValueError(f"{k}는 0~50의 정수입니다.")
            out[k] = int(d[k])
    if d.get("reference") not in (None, ""):
        if d["reference"] != "auto" and d["reference"] not in REFERENCES:
            raise ValueError("reference는 auto, " + ", ".join(REFERENCES) + " 중 하나입니다.")
        out["reference"] = d["reference"]
    prof = str(d.get("profiles") or "").replace(" ", "").upper()
    if prof:
        if not all(re.fullmatch(r"[ACGT]{2,3}>[ACGT]{2,3}", x) for x in prof.split(",")):
            raise ValueError("profiles는 쉼표로 구분한 'GA>AA,TC>TT' 형식입니다.")
        out["profiles"] = prof
    return out


def effective(ds: dict) -> dict:
    return {**defaults_for(ds), **parse_settings(ds.get("mutprofile_json", ""))}


def settings_sha(ds: dict) -> str:
    return hashlib.sha256(json.dumps({**effective(ds), "reference": reference_choice(ds)}, sort_keys=True).encode()).hexdigest()[:16]


# ------------------------------------------------------------------------------------------------ reference + annotation
REFERENCES = {"NC_063383.1": "RefSeq 참조 (clade II 등 기본)", "NC_003310.1": "RefSeq clade I 참조 (Zaire-96-I-16)", "KJ642613.1": "Nextstrain clade-i 빌드 기본 참조 (Congo_8; NC_003310.1과는 다른 유전체)"}


def reference_choice(ds: dict) -> str:
    """The reference accession used: the setting, or - on 'auto' - NC_063383.1 (RefSeq, what Nextclade / Squirrel / snipit use) except for the clade I build,
    which follows the Nextstrain clade-i build's own reference KJ642613.1 (so positions match the Nextstrain tree). NC_003310.1 (RefSeq clade I) can be chosen instead; it is a different genome."""
    r = effective(ds).get("reference", "auto")
    if r in REFERENCES:
        return r
    return "KJ642613.1" if (ds.get("build") or "") == "clade-i" else "NC_063383.1"


def reference_files(ds: dict) -> tuple[Path, Path | None]:
    """The reference FASTA and its gene annotation (GFF3). NC_063383.1 and KJ642613.1 come with the Nextstrain workflow; NC_003310.1 is fetched once from NCBI
    (GenBank record, public) and kept under data/reference/ - the sequence and the genes come from the same record, so they always agree."""
    acc = reference_choice(ds)
    if acc == "NC_003310.1":
        return fetch_genbank_reference(acc)
    base = nextstrain.workflow_dir() / "defaults"
    d = base / "clade-i" if acc == "KJ642613.1" else base
    fa, gff = d / "reference.fasta", d / "genome_annotation.gff3"
    if not fa.exists():
        raise RuntimeError(f"참조 서열({acc})을 찾을 수 없습니다. Nextstrain 설정 카드에서 워크플로를 먼저 받아 주세요.")
    return fa, (gff if gff.exists() else None)


def fetch_genbank_reference(acc: str) -> tuple[Path, Path]:
    from . import config, sources
    d = config.DATA_DIR / "reference"
    fa, gff, gb = d / f"{acc}.fasta", d / f"{acc}.gff3", d / f"{acc}.gb"
    if not (fa.exists() and gff.exists()):
        if not gb.exists():
            d.mkdir(parents=True, exist_ok=True)
            with sources._client() as c:
                r = sources._get(c, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params={"db": "nuccore", "id": acc, "rettype": "gb", "retmode": "text"})
            if r.status_code != 200 or "ORIGIN" not in r.text:
                raise RuntimeError(f"NCBI에서 참조 서열 {acc}를 받지 못했습니다 (HTTP {r.status_code}). 네트워크를 확인하거나 {gb}에 GenBank 파일을 직접 두세요.")
            gb.write_text(r.text, encoding="utf-8")
        seq, gff_text = genbank_to_fasta_gff3(gb.read_text(encoding="utf-8"), acc)
        fa.write_text(f">{acc}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n", encoding="utf-8")
        gff.write_text(gff_text, encoding="utf-8")
    return fa, gff


def _loc_parts(loc: str) -> tuple[list[tuple[int, int]], str]:
    strand = "-" if loc.startswith("complement") else "+"
    parts = [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", loc)]
    if not parts and re.fullmatch(r"(complement\()?\d+\)?", loc):
        n = int(re.search(r"\d+", loc).group())
        parts = [(n, n)]
    return parts, strand


def genbank_to_fasta_gff3(text: str, acc: str) -> tuple[str, str]:
    """The sequence and the CDS features (join / complement included) of one GenBank record as FASTA text and GFF3."""
    seq = "".join(re.findall(r"[acgtn]+", text.split("\nORIGIN", 1)[1].split("\n//")[0])).upper()
    feat = text.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0].split("\n", 1)[1]
    blocks, cur = [], None
    for line in feat.split("\n"):
        if line[5:6].strip():
            cur = [line[5:21].strip(), line[21:].strip(), []]
            blocks.append(cur)
        elif cur is not None:
            t = line[21:].strip()
            if t.startswith("/"):
                cur[2].append(t)
            elif cur[2]:
                cur[2][-1] += " " + t
            else:
                cur[1] += t
    out, n = [f"##gff-version 3", f"##sequence-region {acc} 1 {len(seq)}"], 0
    for key, loc, quals in blocks:
        if key != "CDS":
            continue
        parts, strand = _loc_parts(loc)
        q = {m.group(1): m.group(2) for m in (re.match(r'/(\w+)="?(.*?)"?$', x) for x in quals) if m}
        n += 1
        name = q.get("gene") or q.get("locus_tag") or f"cds{n}"
        attr = f"ID=cds-{n};Name={name};gene={name};product={q.get('product', '').replace(';', '%3B').replace('=', '%3D')};protein_id={q.get('protein_id', '')}"
        for a, b in parts:
            out.append(f"{acc}\tGenBank\tCDS\t{a}\t{b}\t.\t{strand}\t0\t{attr}")
    return seq, "\n".join(out) + "\n"


def read_fasta(path: Path) -> tuple[str, str]:
    name, seq = "", []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if name:
                break
            name = line[1:].split()[0]
        else:
            seq.append(line.strip())
    return name, "".join(seq).upper()


def load_cds(gff: Path | None) -> list[dict]:
    """CDS features: gene name, product, strand, and exon intervals (1-based, inclusive) in coding order."""
    if not gff:
        return []
    by: dict[str, dict] = {}
    for line in open(gff, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        a = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)}
        key = a.get("ID") or a.get("Parent") or f"{f[3]}-{f[4]}"
        c = by.setdefault(key, {"gene": a.get("gene") or a.get("Name") or key, "product": a.get("product", ""), "strand": f[6], "exons": []})
        c["exons"].append((int(f[3]), int(f[4])))
    for c in by.values():
        c["exons"].sort(reverse=(c["strand"] == "-"))
    return list(by.values())


def _coding(c: dict, refseq: str) -> str:
    if c["strand"] == "-":
        return "".join(refseq[a - 1:b].translate(COMP)[::-1] for a, b in c["exons"])
    return "".join(refseq[a - 1:b] for a, b in c["exons"])


def annotate(pos: int, ref: str, alt: str, cds: list[dict], refseq: str) -> dict:
    """Gene / codon / amino acid change of a SNV at 1-based `pos` (reference coordinates). Intergenic -> gene ''. For overlapping genes the first is reported."""
    for c in cds:
        off = 0
        for s, e in c["exons"]:
            if s <= pos <= e:
                idx = off + ((pos - s) if c["strand"] == "+" else (e - pos))      # 0-based position in the coding sequence
                cpos, ci = idx // 3 + 1, idx % 3
                codon = _coding(c, refseq)[cpos * 3 - 3:cpos * 3]
                if len(codon) < 3:
                    break
                a_ = alt if c["strand"] == "+" else alt.translate(COMP)
                new = codon[:ci] + a_ + codon[ci + 1:]
                ra, aa = CODON.get(codon, "?"), CODON.get(new, "?")
                return {"gene": c["gene"], "product": c["product"], "strand": c["strand"], "codon_pos": cpos, "codon_base": ci + 1, "ref_codon": codon, "alt_codon": new,
                        "aa_change": f"{ra}{cpos}{aa}", "effect": "synonymous" if ra == aa else ("stop_gained" if aa == "*" else "nonsynonymous")}
            off += e - s + 1
    return {"gene": "", "product": "", "strand": "", "codon_pos": "", "codon_base": "", "ref_codon": "", "alt_codon": "", "aa_change": "", "effect": "intergenic"}


_CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")


def snvs_from_paf(paf_lines):
    """(query, pos1, ref, alt) for every substitution in minimap2 PAF lines carrying a `cs:Z:` tag. Gaps and N are not reported; a query's SNV seen twice counts once."""
    seen = set()
    for line in paf_lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if not cs:
            continue
        q, rpos = f[0], int(f[7])                          # target start, 0-based
        for t in _CS.findall(cs):
            if t[0] == ":":
                rpos += int(t[1:])
            elif t[0] == "*":
                r, a = t[1].upper(), t[2].upper()
                if r in "ACGT" and a in "ACGT" and (q, rpos + 1) not in seen:
                    seen.add((q, rpos + 1))
                    yield q, rpos + 1, r, a
                rpos += 1
            elif t[0] == "-":
                rpos += len(t) - 1
            elif t[0] == "~":
                rpos += int(re.search(r"\d+", t).group())


# ------------------------------------------------------------------------------------------------ running
def run(ds: dict, cores: int | None = None, force: bool = False, note: str = "", progress=None) -> dict:
    p = progress or (lambda *a: None)
    st = status()
    if not st["ok"]:
        raise RuntimeError(f"mutation-profile / minimap2 환경이 없습니다. 먼저: {st['how']}")
    settings, sha, ver = effective(ds), settings_sha(ds), st["version"]
    members = nextstrain.current_members(ds)
    if not members:
        raise ValueError("데이터셋에 멤버가 없습니다 (서열이 있는 유전체가 있어야 합니다).")
    msha = datasets.members_hash(members)
    if not force:
        for r in nextstrain.list_runs(ds["id"]):
            if r["tool"] == "mutprofile" and r["status"] == "done" and r["members_sha"] == msha and r["params_sha"] == sha and r["workflow_commit"] == ver:
                p(100, f"같은 데이터·설정·버전으로 이미 분석했습니다 ({r['name']})")
                return {**r, "reused": True}
    ref_fa, gff = reference_files(ds)
    refname, refseq = read_fasta(ref_fa)
    cds = load_cds(gff)
    now = datetime.now()
    virus = nextstrain.datasets_virus_name(ds)
    name = nextstrain.slug_name(virus, ds["name"], now.strftime("%Y%m%d-%H%M%S")).replace("-20", "-mutprofile-20", 1)
    taken = {r["name"] for r in nextstrain.list_runs(limit=100000)}
    base, k = name, 1
    while name in taken:
        k += 1
        name = f"{base}-{k}"
    rdir = nextstrain.runs_dir() / export.slug(virus) / name
    (rdir / "mutprofile").mkdir(parents=True, exist_ok=True)
    out, logp = rdir / "mutprofile", rdir / "run.log"
    versions = {"mutation-profile": ver, "minimap2": "", "reference": refname, "annotation": gff.name if gff else ""}
    try:
        versions["minimap2"] = subprocess.run([_bin("minimap2"), "--version"], capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        pass
    con = connect()
    try:
        rid = con.execute("INSERT INTO analysis_run(dataset_id,watch_id,name,status,started,runtime,build,rules_json,params_yaml,params_sha,members_sha,n_public,n_private,workflow_commit,run_dir,log_path,note,versions_json,dataset_name,tool)"
                          " VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                          (ds["id"], ds["watch_id"], name, "running", now.isoformat(timespec="seconds"), "ambient", ds.get("build") or "", json.dumps(ds["rules"], ensure_ascii=False),
                           json.dumps(settings, ensure_ascii=False), sha, msha, sum(1 for m in members if not m.startswith("private:")), sum(1 for m in members if m.startswith("private:")),
                           ver, str(rdir), str(logp), note.strip(), json.dumps(versions), ds["name"], "mutprofile")).lastrowid
        con.commit()
    finally:
        con.close()

    def finish(status_: str, error: str = "", result: str = ""):
        con = connect()
        try:
            con.execute("UPDATE analysis_run SET status=?, finished=?, error=?, auspice_file=? WHERE id=?", (status_, datetime.now().isoformat(timespec="seconds"), error, result, rid))
            audit(con, "mutprofile_run", f"{name}: {status_}" + (f" ({error[:120]})" if error else ""))
            con.commit()
        finally:
            con.close()
    logf = open(logp, "w", encoding="utf-8")
    try:
        p(5, "서열을 모으는 중…")
        fasta = rdir / "input.fasta"
        squirrel.write_input(ds, fasta)
        (rdir / "members.txt").write_text("\n".join(members) + "\n", encoding="utf-8")
        env = tool_env()
        cmd = [_bin("minimap2"), "-x", "asm5", "-c", "--cs", "--secondary=no", "-t", str(cores or 4), str(ref_fa), str(fasta)]
        logf.write("$ " + " ".join(cmd) + "\n")
        p(10, "참조 서열에 정렬하는 중 (minimap2)…")
        paf = rdir / "aln.paf"
        with open(paf, "w") as pf:
            r = subprocess.run(cmd, stdout=pf, stderr=subprocess.PIPE, text=True, env=env)
        logf.write(r.stderr[-3000:])
        if r.returncode != 0:
            raise RuntimeError(f"minimap2가 실패했습니다 (코드 {r.returncode}):\n{r.stderr[-800:]}")
        p(45, "SNP를 읽는 중…")
        snps = list(snvs_from_paf(open(paf, encoding="utf-8")))
        with_snps = len({q for q, *_ in snps})
        n_seq = sum(1 for line in open(fasta, encoding="utf-8") if line.startswith(">"))
        logf.write(f"\nSNP {len(snps):,}건 (서열 {n_seq:,}건 중 SNP가 있는 서열 {with_snps:,}건)\n")
        with open(out / "snps.tsv", "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["ID", "POS", "REF", "ALT"])
            w.writerows(snps)
        if not snps:
            raise RuntimeError("참조 서열과 다른 SNP가 하나도 없습니다 (정렬이 되지 않았을 수 있습니다). 로그를 확인하세요.")
        p(55, "변이 맥락·프로파일을 구하는 중 (mutation-profile)…")
        cmd = [_bin("mutation-profile"), "-f", str(ref_fa), "-m", str(out / "snps.tsv"), "-r", refname, "-b", str(settings["before"]), "-a", str(settings["after"]), "-p", settings["profiles"], "-o", str(out / "profile")]
        logf.write("\n$ " + " ".join(cmd) + "\n")
        r = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=str(rdir))
        logf.write(r.stdout[-3000:] + r.stderr[-3000:])
        if r.returncode != 0 or not (out / "profile.tsv").exists():
            raise RuntimeError(f"mutation-profile이 실패했습니다 (코드 {r.returncode}):\n{(r.stderr or r.stdout)[-800:]}")
        p(80, "유전자·코돈 위치를 붙이는 중…")
        prof = {(row["sample"], int(row["ref_position"])): row for row in csv.DictReader(open(out / "profile.tsv", encoding="utf-8"), delimiter="\t")}
        ann_cache: dict = {}
        genes: dict = defaultdict(Counter)
        sites = set()
        cols = ["sample", "pos", "ref", "alt", "gene", "product", "strand", "aa_change", "effect", "ref_codon", "alt_codon", "codon_pos", "codon_base", "motif_ref", "observed_profile", "profile_of_interest"]
        with open(out / "mutations.tsv", "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(cols)
            for q, pos, r_, a in snps:
                key = (pos, r_, a)
                if key not in ann_cache:
                    ann_cache[key] = annotate(pos, r_, a, cds, refseq)
                an, pr = ann_cache[key], prof.get((q, pos), {})
                row = {"sample": q, "pos": pos, "ref": r_, "alt": a, **an, "motif_ref": pr.get("motif_ref", ""), "observed_profile": pr.get("observed_profile", ""), "profile_of_interest": pr.get("profile_of_interest", "")}
                w.writerow([row[c] for c in cols])
                g = genes[an["gene"] or "(intergenic)"]
                g[an["effect"]] += 1
                g["total"] += 1
                if row["profile_of_interest"] not in ("", "other"):
                    g["poi"] += 1
                sites.add(key)
        with open(out / "gene_summary.tsv", "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["gene", "mutations", "nonsynonymous", "synonymous", "stop_gained", "intergenic", "profile_of_interest"])
            for gname, c in sorted(genes.items(), key=lambda x: -x[1]["total"]):
                w.writerow([gname, c["total"], c["nonsynonymous"], c["synonymous"], c["stop_gained"], c["intergenic"], c["poi"]])
        summary = {"snps": len(snps), "samples_with_snps": with_snps, "samples": n_seq, "distinct_sites": len(sites), "genes_hit": len([g for g in genes if g != "(intergenic)"]),
                   "profile_of_interest": sum(c["poi"] for c in genes.values()), "settings": settings, "reference": refname}
        (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
        logf.write("\n" + json.dumps(summary, ensure_ascii=False) + "\n")
        logf.close()
        fasta.unlink(missing_ok=True)
        paf.unlink(missing_ok=True)
        (rdir / "versions.json").write_text(json.dumps({"tools": versions, "settings": settings}, indent=2), encoding="utf-8")
        finish("done", result="mutprofile/gene_summary.tsv")
        p(100, "분석 완료")
    except BaseException as e:
        if not logf.closed:
            logf.close()
        finish("error", str(e))
        raise
    return nextstrain.get_run(rid)


def result_files(run: dict) -> list[dict]:
    out = Path(run["run_dir"] or "") / "mutprofile"
    if not out.exists():
        return []
    return [{"name": f.name, "path": f"mutprofile/{f.name}", "size": f.stat().st_size} for f in sorted(out.iterdir()) if f.is_file()]


def summary(run: dict) -> dict:
    try:
        return json.loads((Path(run["run_dir"] or "") / "mutprofile" / "summary.json").read_text(encoding="utf-8"))
    except Exception:
        return {}
