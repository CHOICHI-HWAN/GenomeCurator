"""Lineage x mutation heat map (the idea of VIRUS-MVP, github.com/cidgoh/VIRUS-MVP, built on this app's own data) and GVF export.

One row per lineage (or clade), one column per mutation (in genome order, spaced evenly so that every mutation is readable), the colour of a cell is the share of the lineage's genomes
that carry the mutation. Around it: a gene bar, a mutation histogram along the genome, the lineage-defining mutations (from lineage-designation, drawn with a frame), APOBEC3-type marks and
the functional category of the gene. The detail it adds to the snipit-style picture is *per lineage* (frequency instead of single genomes) and *per mutation* (what, where, in whom).

`gvf()` writes the same table as GVF (Genomic Variation Format, the files VIRUS-MVP / nf-ncov-voc read), one file per lineage, so the data can be given to those tools."""
import io
import zipfile
from collections import Counter, defaultdict

from . import insight, lineages, log, mutprofile, outbreak

_log = log.get("heatmap")


def _members(run: dict, group_by: str, max_n: float | None):
    """({genome id: group}, {group: outbreak}) of the run's genomes with a group (lineage / clade) and, with max_n, at most that share of N."""
    ok = mutprofile.good_samples(run, max_n)
    grp, ob = {}, {}
    if group_by == "cluster":                                  # the tree's subclusters (outbreak analysis) as rows
        cl = outbreak.clusters(run, 5, 1, tips=True)
        tip_of = {t: c["id"] for c in cl.get("clusters", []) for t in c["tips"]} if cl.get("available") else {}
        for g in lineages._genomes(run):
            k = tip_of.get(g["id"])
            if k and (ok is None or g["id"] in ok):
                grp[g["id"]] = k
                ob.setdefault(k, g["outbreak"])
        return grp, ob
    for g in lineages._genomes(run):
        k = g["clade"] if group_by == "clade" else g["lineage"]
        if not k or (ok is not None and g["id"] not in ok):
            continue
        grp[g["id"]] = k
        ob.setdefault(k, g["outbreak"])
    return grp, ob


def build(run: dict, group_by: str = "lineage", min_n: int = 3, max_n: float | None = 0.05, min_freq: float = 0.05, min_carriers: int = 2) -> dict:
    grp, ob = _members(run, group_by, max_n)
    rows = insight._mutations(run, max_n)
    ref = mutprofile.summary(run).get("reference", "")
    fa, gff = mutprofile.reference_by_name(ref) if ref else (None, None)
    refname, refseq = mutprofile.read_fasta(fa) if fa else ("", "")
    cds = mutprofile.load_cds(gff)
    size = Counter(grp.values())
    keep = {g for g, n in size.items() if n >= min_n}
    carriers: dict[tuple, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    info: dict[tuple, dict] = {}
    for r in rows:
        g = grp.get(r["sample"])
        if g not in keep:
            continue
        k = (int(r["pos"]), r["alt"])
        carriers[k][g].add(r["sample"])
        info.setdefault(k, r)
    sites = []
    for k, per in carriers.items():
        if sum(len(v) for v in per.values()) < min_carriers:
            continue
        if max(len(v) / size[g] for g, v in per.items()) < min_freq:
            continue
        sites.append(k)
    sites.sort()
    # lineage-defining mutations in the run's coordinates, per group
    defset: dict[str, dict] = {}
    to_run = None
    have_defs = False
    try:
        if lineages.status()["ok"] and fa:
            to_run = lineages.defining_in_run_coordinates(run, fa)
            have_defs = True
    except Exception:
        _log.warning("lineage definitions unavailable for the heat map", exc_info=True)
    if have_defs:
        for g in keep:
            if group_by == "lineage" and ob.get(g):
                defset[g] = {((to_run.get(p) if to_run is not None else p), n): by for p, n, by in lineages.defining_snps(ob[g], g)}
    defined_by: dict[tuple, str] = {}
    for g, d in defset.items():
        for key, by in d.items():
            if key[0]:
                defined_by.setdefault(key, by)
    cat = outbreak.categories()
    cols = []
    for j, k in enumerate(sites):
        r = info[k]
        c = cat.get(r["gene"], {}) if r["gene"] else {}
        cols.append({"label": f"{r['ref']}{k[0]}{k[1]}", "pos": k[0], "ref": r["ref"], "alt": k[1], "gene": r["gene"], "aa": r["aa_change"], "effect": r["effect"], "apobec": insight._is_apo(r),
                     "n": sum(len(v) for v in carriers[k].values()), "defined_by": defined_by.get(k, ""), "function": c.get("function", ""), "category": c.get("category", ""),
                     "profile": r.get("profile_of_interest", "")})
    meta = {g["id"]: g for g in lineages._genomes(run)}
    order = []
    for g in keep:
        ds = sorted(d for d in (lineages.parse_date(meta[i]["date"])[3] for i, gg in grp.items() if gg == g and i in meta) if d)
        order.append((ds[0] if ds else "9999", g))
    order.sort()
    rows_out, cells = [], []
    for ri, (_f, g) in enumerate(order):
        n = size[g]
        gs = [meta[i] for i, gg in grp.items() if gg == g and i in meta]
        ds = sorted(d for d in (lineages.parse_date(x["date"])[3] for x in gs) if d)
        rows_out.append({"id": g, "n": n, "outbreak": ob.get(g, ""), "first": ds[0] if ds else "", "last": ds[-1] if ds else "", "countries": dict(Counter(x["country"] for x in gs if x["country"]).most_common(4)),
                         "designated": bool(defset.get(g)), "n_defining": len(defset.get(g, {}))})
        for j, k in enumerate(sites):
            c = len(carriers[k].get(g, ()))
            if c:
                cells.append([ri, j, round(100 * c / n, 1), c, 1 if k in defset.get(g, {}) else 0])
    L = len(refseq)
    hist = Counter(k[0] // 1000 for k in sites)
    genes = [{"gene": c["gene"], "start": min(a for a, _ in c["exons"]), "end": max(b for _, b in c["exons"]), "strand": c["strand"], "category": cat.get(c["gene"], {}).get("category", "")} for c in cds]
    win = mutprofile.run_itr(run)
    return {"reference": refname or ref, "length": L, "group_by": group_by, "rows": rows_out, "columns": cols, "cells": cells, "genes": genes, "histogram": [[b * 1000, n] for b, n in sorted(hist.items())],
            "itr": ([[1, win[0] - 1], [win[1] + 1, L]] if win else None), "definitions": have_defs, "max_n": max_n, "min_n": min_n, "min_freq": min_freq,
            "n_genomes": len(grp), "n_dropped_groups": sum(1 for n in size.values() if n < min_n)}


def gvf(run: dict, group: str = "", group_by: str = "lineage", max_n: float | None = 0.05) -> bytes | str:
    """GVF text of one lineage (or clade), or a zip of one GVF per lineage with at least 3 genomes when `group` is empty. Columns after the standard nine: the attributes the VIRUS-MVP /
    nf-ncov-voc GVFs carry that this app knows (mutation description, gene, protein change, reference / alternate counts, depth, sample size) and `clade_defining`."""
    h = build(run, group_by, 1, max_n, 0.0, 1)
    rows, cols = h["rows"], h["columns"]
    by_row = defaultdict(dict)
    for ri, cj, _f, c, d in h["cells"]:
        by_row[ri][cj] = (c, d)
    tax = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10244"

    def one(ri: int) -> str:
        r = rows[ri]
        out = ["##gff-version 3", "##gvf-version 1.10", f"##sample-description sample_desc=Clinical;sample_group={r['id']};", f"##species {tax}",
               "#seqid\t#source\t#type\t#start\t#end\t#score\t#strand\t#phase\t#attributes"]
        for i, (cj, (c, d)) in enumerate(sorted(by_row[ri].items())):
            m = cols[cj]
            desc = m["aa"] if m["aa"] and m["gene"] else f"g.{m['ref']}{m['pos']}{m['alt']}"
            attrs = {"ID": f"ID_{i}", "original_mutation_description": desc, "gene": m["gene"] or "intergenic", "product": m["function"] or "n/a", "ro": r["n"] - c, "ao": c, "dp": r["n"],
                     "sample_size": r["n"], "Reference_seq": m["ref"], "Variant_seq": m["alt"], "clade_defining": "True" if d else "False", "apobec3_type": "True" if m["apobec"] else "False"}
            out.append(f"{h['reference'].split()[0]}\t.\tSNV\t{m['pos']}\t{m['pos']}\t.\t+\t.\t" + ";".join(f"{k}={str(v).replace(';', ',')}" for k, v in attrs.items()))
        return "\n".join(out) + "\n"
    if group:
        ri = next((i for i, r in enumerate(rows) if r["id"] == group), None)
        if ri is None:
            raise ValueError(f"{group}을(를) 찾을 수 없습니다.")
        return one(ri)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for ri, r in enumerate(rows):
            if r["n"] >= 3:
                z.writestr(f"{r['id']}_annotated.gvf", one(ri))
    return buf.getvalue()
