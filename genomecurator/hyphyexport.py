"""Input for site-wise selection tests (HyPhy FEL / MEME / FUBAR): a codon alignment per gene and the tree (evidence E10 of docs/ADAPTIVE.md).

The app does not run these tests - validated programs exist - it only prepares what they need from the run's own data:
* **genes/<gene>.fasta**: the coding sequence of every genome in the gene, cut from Squirrel's alignment (all rows are in the alignment reference's coordinates; the reference's columns are mapped to the run's
  reference, whose annotation gives the exons). A gene is skipped (and listed in the README) when a codon position of the run's reference has no aligned column. Deleted bases are `-`, masked bases `N`.
  Genomes that carry an in-frame stop codon are kept (HyPhy needs them cleaned: see its pre-msa step), deletions that are not a multiple of three break the frame and are reported per gene.
* **tree.nwk**: the dataset's Nextstrain tree as Newick, branch length = substitutions per site (the change of `div` divided by the genome length), restricted to the genomes of the alignment.
Only public genomes are written (a private genome's id must not leave the app in a file that is meant to be shared)."""
import io
import zipfile

from . import aligned, mutprofile, nextstrain, outbreak, squirrel

_COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")


def gene_positions(c: dict) -> list[tuple[int, str]]:
    """The reference positions of a CDS in coding order, with the strand ('+' / '-') used to read the base."""
    out = []
    for a, b in c["exons"]:
        rng = range(a, b + 1) if c["strand"] == "+" else range(b, a - 1, -1)
        out += [(p, c["strand"]) for p in rng]
    return out


def codon_alignment(rows: dict[str, str], inv: dict[int, int], cds: dict) -> tuple[dict[str, str] | None, str]:
    """{genome: coding sequence} for one CDS from alignment rows, or (None, reason). `inv` maps the run reference's positions to alignment columns (1-based)."""
    pos = gene_positions(cds)
    if len(pos) % 3:
        return None, "CDS 길이가 3의 배수가 아님"
    cols = []
    for p, _st in pos:
        q = inv.get(p)
        if q is None:
            return None, "정렬 열이 없는 코돈 위치가 있음"
        cols.append(q - 1)
    minus = cds["strand"] == "-"
    out = {}
    for name, seq in rows.items():
        if cols and max(cols) >= len(seq):
            return None, "정렬 행이 참조보다 짧음"
        s = "".join(seq[i] for i in cols).upper()
        out[name] = s.translate(_COMP) if minus else s
    return out, ""


def newick(tree: dict, keep: set[str], length: int) -> str:
    """Newick of the tree restricted to the tips in `keep`; branch length = (div - parent div) / genome length."""
    import sys
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))

    def walk(n, pdiv):
        div = (n.get("node_attrs", {}).get("div") or pdiv)
        bl = max(0.0, (div - pdiv)) / max(1, length)
        kids = [k for k in (walk(c, div) for c in (n.get("children") or [])) if k]
        if not n.get("children"):
            return f"{n['name']}:{bl:.8f}" if n["name"] in keep else ""
        if not kids:
            return ""
        return f"({','.join(kids)}):{bl:.8f}"
    return walk(tree, 0.0) + ";"


def build(run: dict, genes: list[str] | None = None, max_genes: int = 60) -> bytes:
    """A zip with README.txt, tree.nwk and genes/<gene>.fasta for the requested genes (default: none = all, capped at `max_genes`)."""
    ref = mutprofile.summary(run).get("reference", "")
    fa, gff = mutprofile.reference_by_name(ref)
    cds_all = {c["gene"]: c for c in mutprofile.load_cds(gff)}
    if genes:
        miss = [g for g in genes if g not in cds_all]
        if miss:
            raise ValueError(f"참조 주석에 없는 유전자: {', '.join(miss)}")
        want = genes
    else:
        want = sorted(cds_all)[:max_genes]
    sm = mutprofile.summary(run)
    from . import datasets
    ds = datasets.get(run["dataset_id"])
    arun = next((r for r in nextstrain.list_runs(run["dataset_id"]) if r["tool"] == "squirrel" and r["status"] == "done" and squirrel.alignment_path(r)), None)
    aln = squirrel.alignment_path(arun) if arun else None
    aref = squirrel.alignment_reference(ds)
    if not aln or not aref:
        raise ValueError("Squirrel 정렬이 필요합니다(이 데이터셋의 Squirrel 분석을 먼저 실행하세요).")
    afa = squirrel.reference_fasta(aref)
    members = {m for m in nextstrain.run_members(run) if not str(m).startswith("private:")}
    rows = {n: s for n, s in aligned.read_alignment(aln) if n in members}
    cmap = aligned._liftover(afa, fa) if str(afa) != str(fa) else {i: i for i in range(1, len(mutprofile.read_fasta(fa)[1]) + 1)}
    inv = {t: q for q, t in cmap.items()}
    tree, _root, _ns = outbreak._tree(run)
    notes, written = [], 0
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for g in want:
            seqs, why = codon_alignment(rows, inv, cds_all[g])
            if seqs is None:
                notes.append(f"{g}: 건너뜀 - {why}")
                continue
            ln = len(next(iter(seqs.values()))) if seqs else 0
            bad = sum(1 for s in seqs.values() if any(s[i:i + 3] in ("TAA", "TAG", "TGA") for i in range(0, max(0, ln - 3), 3)))
            gap = sum(1 for s in seqs.values() if s.count("-") % 3)
            notes.append(f"{g}: 서열 {len(seqs)}건, 길이 {ln}bp" + (f", 안쪽 종결 코돈이 있는 서열 {bad}건" if bad else "") + (f", 3의 배수가 아닌 결실이 있는 서열 {gap}건" if gap else ""))
            z.writestr(f"genes/{g}.fasta", "".join(f">{n}\n{s}\n" for n, s in seqs.items()))
            written += 1
        if tree:
            z.writestr("tree.nwk", newick(tree, set(rows), len(mutprofile.read_fasta(fa)[1])) + "\n")
        z.writestr("README.txt", "HyPhy 입력 (Virus Genome Metadata Curator)\n\n실행: " + run["name"] + f"\n참조: {ref} · 정렬: Squirrel({aref})\n공개 유전체 {len(rows)}건(비공개는 뺌)\n\n"
                   "genes/<유전자>.fasta = 코돈 정렬(코딩 방향), tree.nwk = Nextstrain 트리(가지 길이 = 자리당 치환 수).\n"
                   "HyPhy(FEL / MEME / FUBAR)는 이 앱이 실행하지 않습니다. 종결 코돈·3의 배수가 아닌 결실은 HyPhy의 pre-msa 단계로 정리하세요.\n"
                   "이 파일들에는 공개 유전체의 accession이 들어 있습니다. 공유 전에 확인하세요.\n\n" + "\n".join(notes) + "\n")
    return buf.getvalue()
