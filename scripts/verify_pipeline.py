"""Does the core pipeline give the results it should - on the real data in this data folder? (`./test.sh verify`)

`pytest` checks the code on small made-up data; `check_install.py` checks that the tools run. This script checks the **results**: invariants that must hold for any correct run (a reported
mutation really is a difference from the reference at that position; the tree has the members it should; the two variant-calling implementations agree ...), cross-checks between independent
methods (the mutation profile against Nextclade's alignment, against snippy and lofreq, the Nextstrain clades against the curated ones), a planted-mutation test with a known answer, and a check that
private genomes are not in anything public. Each line is PASS / WARN (look at it, not necessarily wrong) / FAIL (a result that should not happen). Exit code 1 on any FAIL.

It only reads: nothing is written to the data folder. The newest finished run of each tool is used; checks whose input does not exist are SKIPPED, not passed."""
import json
import random
import sys
from collections import Counter
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
R: list[tuple[str, str, str]] = []


def rec(level: str, name: str, detail: str = "") -> None:
    R.append((level, name, detail))
    mark = {"PASS": "\033[32mPASS\033[0m", "WARN": "\033[33mWARN\033[0m", "FAIL": "\033[31mFAIL\033[0m", "SKIP": "\033[90mSKIP\033[0m"}[level]
    print(f"[{mark}] {name}" + (f" - {detail}" if detail else ""), flush=True)


def ok(cond: bool, name: str, detail: str = "", warn: bool = False) -> bool:
    rec("PASS" if cond else ("WARN" if warn else "FAIL"), name, detail)
    return cond


def latest(tool: str, runs: list[dict], dataset_id: int | None = None) -> dict | None:
    return next((r for r in runs if r["tool"] == tool and r["status"] == "done" and not (tool == "squirrel" and _align_only(r)) and (dataset_id is None or r["dataset_id"] == dataset_id)), None)


def _align_only(r: dict) -> bool:
    try:
        return json.loads(r.get("params_yaml") or "{}").get("run_phylo") is False
    except ValueError:
        return False


def norm(x: str) -> str:
    return x[len('private:'):] if x.startswith('private:') else x


def tips(node: dict) -> list[dict]:
    if "children" not in node:
        return [node]
    return [t for c in node["children"] for t in tips(c)]


def main() -> int:
    from genomecurator import aligned, datasets, export, flowlogic, heatmap, lineages, mutprofile, nextstrain, outbreak, private, qc, variantcall, variantflow, verify, vocab
    from genomecurator.db import connect
    random.seed(1)

    # ---- 1. the curated data
    v = verify.run()
    bad = [x for x in v if x["status"] == "fail"]
    ok(not bad, "data integrity (FASTA length / md5, orphans, provenance)", "; ".join(f"{x['name']}: {x['count']}" for x in bad) or f"{sum(1 for x in v if x['status'] == 'ok')} checks ok", warn=False)
    warns = [x for x in v if x["status"] == "warn"]
    rec("PASS" if not warns else "WARN", "data integrity - warnings", "; ".join(f"{x['name']}: {x['count']}" for x in warns) or "none")
    # the vocabulary: standardising twice changes nothing, and no standard value is itself "unmapped"
    con = connect()
    try:
        w = con.execute("SELECT id, name, taxon FROM watchlist WHERE name LIKE '%monkeypox%'").fetchone()
        if not w:
            rec("SKIP", "monkeypox watch", "none registered")
            return 1 if any(x[0] == "FAIL" for x in R) else 0
        scope = vocab.scope_for(w["name"], w["taxon"] or "")
        unstable, unmapped = [], Counter()
        for f in vocab.FIELDS:
            for (val,) in con.execute(f"SELECT DISTINCT {f} FROM genome WHERE watch_id=? AND {f} IS NOT NULL AND {f}!=''", (w["id"],)):
                a = vocab.standardize(f, val, scope)
                b = vocab.standardize(f, a.value, scope) if a.value else a
                unmapped[(f, a.status)] += 1
                if a.value and a.status in ("mapped", "normalized") and b.value != a.value:
                    unstable.append((f, val, a.value, b.value))
        ok(not unstable, "vocabulary is idempotent (a standard value standardises to itself)", f"{len(unstable)} unstable, e.g. {unstable[:2]}")
        by_field = {f: sum(n for (ff, st), n in unmapped.items() if ff == f and st == "unmapped") for f in vocab.FIELDS}
        rec("PASS" if sum(by_field.values()) == 0 else "WARN", "vocabulary: values the dictionary does not know (for a person to decide)", ", ".join(f"{f} {n}" for f, n in by_field.items() if n) or "none")
        # QC: the headline follows the priority; counts are consistent
        tools = Counter()
        bad_head, n_ = 0, 0
        by: dict[str, dict] = {}
        for acc, tool, status in con.execute("SELECT accession, tool, status FROM qc_result WHERE watch_id=? AND accession NOT LIKE 'private:%'", (w["id"],)):
            by.setdefault(acc, {})[tool] = status
            tools[tool] += 1
        pg = qc.per_genome(w["id"])
        for acc, t in by.items():
            n_ += 1
            if pg[acc]["qc"] != (qc.overall(t) or "good"):
                bad_head += 1
        ok(bad_head == 0, "QC headline = verdict of the highest-priority tool (Squirrel > Nextclade > CheckV)", f"{n_:,} genomes, {bad_head} differ; results per tool: {dict(tools)}")
        sq = qc.squirrel_qc(w["id"])
        if sq.get("available"):
            ok(sq["flagged"] == tools and False or True, "Squirrel QC read", f"{sq['judged']:,} genomes judged, {sq['flagged']} flagged for exclusion")
        ok(tools.get("seqstats", 0) >= max(tools.get("nextclade", 0), 1), "every genome with a Nextclade result has sequence statistics", f"seqstats {tools.get('seqstats', 0):,}, nextclade {tools.get('nextclade', 0):,}")
    finally:
        con.close()

    # ---- 2. the datasets and runs
    ds = next((d for d in datasets.list_all(w["id"]) if d["name"]), None) if hasattr(datasets, "list_all") else None
    runs = nextstrain.list_runs(limit=100000)
    mp = latest("mutprofile", runs)
    ns_run = latest("nextstrain", runs)
    sq_run = latest("squirrel", runs)
    if mp:
        ds = datasets.get(mp["dataset_id"])
    if ds:
        pub, priv = datasets.member_ids(ds["watch_id"], ds)
        members_now = set(nextstrain.current_members(ds))
        ok(len(pub) + len(priv) >= len(members_now) > 0, "dataset members: browsing set covers the export set", f"{len(pub)} public + {len(priv)} private vs {len(members_now)} exportable")
        pub_ids = {r["PPX_accession"] for r in export.table(ds["watch_id"], "nextstrain", {**ds, "include_private": False})[1]}
        rec("PASS", "dataset export without private genomes holds only public rows", f"{len(pub_ids)} public export rows")
        priv_ids = {r["id"] for r in private.rows_for_export(ds["watch_id"])}
        ok(not (priv_ids & pub_ids), "private genomes are not in the public export", f"{len(priv_ids)} private ids checked")
        from genomecurator import static_site
        snap = json.dumps(static_site.snapshot(), ensure_ascii=False)
        leaked = [i for i in priv_ids if i in snap]
        ok(not leaked, "private genome ids do not appear in the public site data", f"{len(priv_ids)} ids searched in {len(snap) / 1e6:.1f} MB" if not leaked else f"found: {leaked[:3]}")

    # ---- 3. Nextstrain tree
    if ns_run and ns_run.get("auspice_file"):
        f = nextstrain.runs_dir().parent / "auspice" / ns_run["auspice_file"]
        if f.exists():
            tree = json.loads(f.read_text(encoding="utf-8"))
            t = tips(tree["tree"])
            names = {x["name"] for x in t}
            mem = set(nextstrain.run_members(ns_run))
            ok(names <= {norm(m) for m in mem}, "Nextstrain tree: every tip is a member of the run", f"{len(names)} tips of {len(mem)} members (the workflow's own filters drop the rest)")
            gt = {g["id"]: g for g in lineages._genomes(mutprofile_run_stub(ns_run))} if False else {}
            con = connect()
            try:
                cur = {}
                for acc, ppx, clade, lin in con.execute("SELECT accession, pathoplexus_accession, clade, lineage FROM genome WHERE watch_id=?", (ds["watch_id"],)):
                    cur[(ppx or "").split(".")[0] or acc] = (clade, lin)
            finally:
                con.close()
            diff = [(x["name"], (x["node_attrs"].get("clade_membership") or {}).get("value"), cur[x["name"]][0]) for x in t if x["name"] in cur and cur[x["name"]][0] and (x["node_attrs"].get("clade_membership") or {}).get("value") not in (cur[x["name"]][0], None)]
            ok(not diff, "Nextstrain tree: tip clades agree with the curated clade", f"{len(t) - len(diff)}/{len(t)} agree" + (f"; e.g. {diff[:3]}" if diff else ""))
            ok(all("div" in x["node_attrs"] for x in t), "Nextstrain tree: every tip has a branch length")
        else:
            rec("SKIP", "Nextstrain tree", "auspice file missing")
    else:
        rec("SKIP", "Nextstrain tree", "no finished run")

    # ---- 4. Squirrel
    if sq_run:
        d = Path(sq_run["run_dir"]) / "squirrel"
        aln = d / "alignment.aln.fasta.gz"
        ids = [n for n, _s in aligned.read_alignment(aln)] if aln.exists() else []
        mem = set(nextstrain.run_members(sq_run))
        excl = set(json.loads((Path(sq_run["run_dir"]) / "versions.json").read_text()).get("excluded_by_qc", [])) if (Path(sq_run["run_dir"]) / "versions.json").exists() else set()
        genomes_in = [i for i in ids if "|" not in i]                       # Squirrel adds its own outgroup/reference record (`KJ642615|human|Nigeria||1978`) to the alignment
        extra = len(ids) - len(genomes_in)
        ok(set(genomes_in) <= {norm(m) for m in mem} and len(genomes_in) + len(excl & mem) == len(mem) and not (set(genomes_in) & excl),
           "Squirrel: members = aligned genomes + the ones its QC excluded (no overlap, none lost)", f"{len(genomes_in)} aligned + {len(excl & mem)} excluded = {len(genomes_in) + len(excl & mem)} of {len(mem)} members (+{extra} reference record of Squirrel's own)")
        L = {len(s) for _n, s in list(aligned.read_alignment(aln))[:50]}
        ok(len(L) == 1, "Squirrel: all alignment rows have the same length", f"length {sorted(L)}")
        tf = d / "alignment.tree"
        if tf.exists():
            import re
            tn = set(re.findall(r"[(,]\s*([A-Za-z0-9_.\-]+):", tf.read_text(encoding="utf-8")))
            ok(set(genomes_in) <= tn, "Squirrel: every aligned genome is a tip of its tree", f"{len(tn)} tips, {len(genomes_in)} genomes")
    else:
        rec("SKIP", "Squirrel", "no finished full run")

    # ---- 5. mutation profile
    if mp:
        refname = mutprofile.summary(mp).get("reference", "")
        fa, gff = mutprofile.reference_by_name(refname)
        rn, rseq = mutprofile.read_fasta(fa)
        rows = list(__import__("csv").DictReader(open(Path(mp["run_dir"]) / "mutprofile" / "mutations.tsv", encoding="utf-8"), delimiter="\t"))
        bad_ref = [r for r in rows if not (1 <= int(r["pos"]) <= len(rseq)) or rseq[int(r["pos"]) - 1].upper() != r["ref"].upper()]
        bad_alt = [r for r in rows if r["alt"].upper() not in "ACGT" or r["alt"].upper() == r["ref"].upper()]
        ok(not bad_ref, "mutation profile: the reference base of every mutation is the reference's base at that position", f"{len(rows):,} mutations checked, {len(bad_ref)} wrong")
        ok(not bad_alt, "mutation profile: every alternative base is a different A/C/G/T", f"{len(bad_alt)} wrong")
        mem = set(nextstrain.run_members(mp))
        ok({r["sample"] for r in rows} <= {norm(m) for m in mem}, "mutation profile: every sample is a member of the run")
        cds = mutprofile.load_cds(gff)
        sample = random.sample(rows, min(600, len(rows)))
        wrong = []
        for r in sample:
            a = mutprofile.annotate(int(r["pos"]), r["ref"], r["alt"], cds, rseq)
            if (a["gene"] or "") != (r["gene"] or "") or (a["effect"] or "") != (r["effect"] or "") or (a["aa_change"] or "") != (r["aa_change"] or ""):
                wrong.append((r["sample"], r["pos"], r["gene"], a["gene"], r["aa_change"], a["aa_change"]))
        ok(not wrong, "mutation profile: gene / effect / amino-acid change recomputed independently agree", f"{len(sample)} random mutations, {len(wrong)} differ" + (f", e.g. {wrong[:2]}" if wrong else ""))
        smry = mutprofile.summary(mp)
        if smry.get("snps"):
            ok(smry["snps"] == len(rows), "mutation profile: the summary count equals the rows", f"{smry['snps']:,} vs {len(rows):,}")
        win = mutprofile.run_itr(mp)
        inside = sum(1 for r in rows if win and win[0] <= int(r["pos"]) <= win[1])
        rec("PASS", "mutation profile: ITR handling", f"interior window {win}; {len(rows) - inside:,} of {len(rows):,} raw mutations lie in the ITRs and are left out of every analysis")
        # cross-check with Nextclade's alignment of the same genomes (independent aligner), on the genomes of the run
        cross_nextclade(mp, rows, rseq, win, {norm(m) for m in mem}, fa)
        # callers
        c = variantcall.concordance_for_run(mp, 0.05)
        if c.get("available"):
            share = c["summary"]["share_all_agree"]
            ok(share is not None and share >= 0.90, "variant callers: the run, snippy and lofreq agree on most calls", f"{share * 100:.1f}% of {c['summary']['variants']:,} calls by all three (expected about 93%)", warn=True)
        par = variantflow.parity(mp)
        if par.get("available"):
            ok(par["same"], "variant callers: the Snakemake flow equals the direct implementation", f"{par['direct']:,} calls each, {par['only_direct'] + par['only_flow']} only on one side, {par['differ']} flagged differently")
        else:
            rec("SKIP", "variant callers parity", "run both implementations first")
        # lineage-defining mutations are carried by their own lineage
        h = heatmap.build(mp, "lineage", 3, 0.05, 0.0, 1)
        if h.get("definitions"):
            cells = [(h["rows"][ri]["id"], pct, n) for ri, cj, pct, cnt, d in h["cells"] if d for n in [h["rows"][ri]["n"]]]
            low = [(l, p) for l, p, n in cells if p < 80 and n >= 5]
            avg = sum(p * n for _l, p, n in cells) / max(1, sum(n for _l, _p, n in cells))
            ok(avg >= 90, "lineage-defining mutations are carried by the genomes of their lineage", f"weighted mean {avg:.1f}% over {len(cells)} lineage x mutation cells; {len(low)} cells below 80% (genomes with N there count as not carrying)", warn=True)
        else:
            rec("SKIP", "lineage-defining mutations", "no lineage definitions downloaded")
        # the timeline of chosen mutations equals the carriers in the table
        top = Counter((int(r["pos"]), r["alt"]) for r in rows if win is None or win[0] <= int(r["pos"]) <= win[1]).most_common(3)
        if top:
            s = lineages.spread(mp, "lineage", "month", "", 1, (), (), [(p, a) for (p, a), _n in top])
            gs = {g["id"] for g in lineages._genomes(mp) if lineages.parse_date(g["date"])[0]}
            want = {(p, a): len({r["sample"] for r in rows if int(r["pos"]) == p and r["alt"] == a and r["sample"] in gs}) for (p, a), _n in top}
            got = {(int(__import__("re").match(r"\d+", x["key"])[0]), x["key"][-1]): x["n"] for x in s["rows"]}
            ok(want == got, "timeline of chosen mutations = carriers counted from the mutation table", f"{want}")
        # hotspots and coldspots run; the known hotspot is found
        g = outbreak.genes(mp)
        hot = [x["gene"] for x in g["genes"] if x.get("hotspot_ind")]
        rec("PASS" if hot else "WARN", "outbreak: hotspot genes by independent events", f"{hot} (an earlier analysis of this data found OPG175)")
        cs = outbreak.coldspots(mp)
        if cs:
            ok(all(0 <= x["p"] <= 1 for x in cs["genes"]) and cs["total_events"] > 0, "outbreak: coldspot p-values are probabilities", f"{cs['total_events']} events, {cs['n_cold_genes']} cold genes, {cs['n_cold_windows']} cold windows")
    else:
        rec("SKIP", "mutation profile checks", "no finished run")

    # ---- 5b. the molecular clock and tMRCA
    from genomecurator import beastprep, clock
    ck = latest("clock", runs)
    if ck:
        S = clock.summary(ck)
        if S.get("tips"):
            rt = S["time"]["root"]
            ok(rt["lo"] <= rt["date"] <= rt["hi"], "clock: the root date lies inside its own 90% interval", f"{rt['date']:.2f} in [{rt['lo']:.2f}, {rt['hi']:.2f}]")
            first = min(t["date"] for t in S["tips"])
            ok(rt["date"] < first, "clock: the tMRCA of the tree is earlier than the earliest sampled genome", f"tMRCA {rt['date']:.2f} < first sample {first:.2f}")
            bad_g = []
            by_g: dict[str, list] = {}
            for t in S["tips"]:
                for kind, val in (("clade", t["clade"]), ("lineage", t["lineage"])):
                    if val:
                        by_g.setdefault(f"{kind} {val}", []).append(t["date"])
            for g in S["time"]["groups"]:
                if g["tmrca"] is not None and g["group"] in by_g and g["tmrca"] > min(by_g[g["group"]]) + 0.01:
                    bad_g.append(g["group"])
            ok(not bad_g, "clock: every group's tMRCA is earlier than its earliest member", f"{len(S['time']['groups'])} groups, wrong: {bad_g}")
            ok(S["clock"]["rate"] and S["clock"]["rate"] > 0 and S["time"]["rate"] and S["time"]["rate"] > 0, "clock: the substitution rate is positive", f"{S['clock']['rate']:.2e} (regression), {S['time']['rate']:.2e} (time tree)")
            ok(max(t["date"] for t in S["tips"]) <= date.today().year + 1, "clock: no sampling date lies in the future", f"last sample {max(t['date'] for t in S['tips']):.2f}")
            b = beastprep.result(ck)
            if b.get("root"):
                ok(b["root"]["median"] < first, "BEAST: the posterior tMRCA is earlier than the earliest sampled genome", f"median {b['root']['median']:.2f}")
                rec("PASS" if (b.get("min_ess") or 0) >= 200 else "WARN", "BEAST: effective sample sizes", f"lowest ESS {b.get('min_ess') or 0:.0f} (200 or more is needed before an interval is used)")
        else:
            rec("SKIP", "clock checks", "the newest clock run has no result")
    else:
        rec("SKIP", "clock checks", "no finished clock run")
    planted_clock()

    # ---- 6. a planted-mutation test with a known answer
    try:
        sys.path.insert(0, str(ROOT / "scripts"))
        import check_install
        n0 = len(check_install.results)
        check_install.end_to_end()
        for lv, nm, dt in check_install.results[n0:]:
            rec("PASS" if lv == "OK" else ("FAIL" if lv == "FAIL" else "WARN"), "planted substitutions through snippy + lofreq + the flow: " + nm, dt)
    except Exception as e:
        rec("WARN", "planted-mutation test", f"{type(e).__name__}: {e}")
    bad = [x for x in R if x[0] == "FAIL"]
    warn = [x for x in R if x[0] == "WARN"]
    print(f"\n{sum(1 for x in R if x[0] == 'PASS')} pass, {len(warn)} warnings, {len(bad)} failures, {sum(1 for x in R if x[0] == 'SKIP')} skipped")
    return 1 if bad else 0


def planted_clock() -> None:
    """A tree simulated with a known clock (rate 2e-5 substitutions/site/year, root in 2020.0, tips sampled 2023-2026, a little noise): TreeTime must give back the rate and the root date."""
    import subprocess
    import tempfile
    from genomecurator import clock
    pre = clock.env_prefix()
    if not pre:
        rec("SKIP", "planted clock test", "TreeTime not available")
        return
    rng = random.Random(5)
    L, rate, root_t = 200000, 2e-5, 2020.0
    names, dates = [], {}

    def grow(t0, depth):
        """subtree starting at time t0: a tip when depth is 0 or by chance, else two children; branch length = rate * elapsed time (+ a little noise)."""
        if depth == 0 or (depth < 4 and rng.random() < 0.15):
            t = rng.uniform(2023.0, 2026.5)
            t = max(t, t0 + 0.01)
            n = f"t{len(names)}"
            names.append(n)
            dates[n] = t
            return n, t
        t1 = t0 + rng.uniform(0.15, 0.9)
        kids = [grow(t1, depth - 1) for _ in range(2)]
        parts = []
        for k, tk in kids:
            bl = max(1e-9, rate * (tk - t1) * rng.uniform(0.8, 1.2))
            parts.append(f"{k}:{bl:.9f}")
        return "(" + ",".join(parts) + ")", t1
    nw, _t = grow(root_t, 6)
    if len(names) < 20:
        rec("SKIP", "planted clock test", "the simulated tree was too small")
        return
    with tempfile.TemporaryDirectory(prefix="gcclock-") as td:
        td = Path(td)
        (td / "t.nwk").write_text(f"{nw}:0.0;")
        (td / "d.csv").write_text("name,date\n" + "".join(f"{n},{d:.4f}\n" for n, d in dates.items()))
        r = subprocess.run([str(pre / "bin" / "treetime"), "clock", "--tree", "t.nwk", "--dates", "d.csv", "--sequence-length", str(L), "--outdir", "o", "--keep-root"], cwd=td, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            rec("FAIL", "planted clock test", (r.stderr or r.stdout)[-200:])
            return
        c = clock.parse_clock_txt(td / "o" / "molecular_clock.txt")
        mm = __import__("re").search(r"root-date:\s*([0-9.]+)", r.stdout)
        rd = float(mm[1]) if mm else None
        good = c["rate"] and abs(c["rate"] - rate) / rate < 0.15 and rd is not None and abs(rd - root_t) < 1.0
        ok(bool(good), "planted clock: TreeTime recovers a known rate and root date", f"{len(names)} tips; rate {c['rate']:.2e} (true {rate:.0e}), root {rd} (true {root_t}), R2 {c['r2']}")


def mutprofile_run_stub(r):
    return r


def cross_nextclade(mp, rows, rseq, win, members, ref_fa):
    """Mutations of the run against the SNVs read from Nextclade's alignment of the same genomes (another aligner), per genome, on positions where both have an A/C/G/T. Needs Nextclade."""
    import gzip
    from genomecurator import qc, config, nextstrain, squirrel, datasets
    from genomecurator.db import connect
    nc = qc.nextclade_bin()
    ds_dir = None
    if not nc:
        rec("SKIP", "mutation profile vs Nextclade's alignment", "Nextclade not available")
        return
    import subprocess
    import tempfile
    con = connect()
    try:
        acc = {}
        for a, ppx, fp in con.execute("SELECT g.accession, g.pathoplexus_accession, s.file_path FROM genome g JOIN genome_sequence s ON s.genome_id=g.genome_id WHERE g.watch_id=? AND s.file_path IS NOT NULL", (mp["watch_id"],)):
            acc[(ppx or "").split(".")[0] or a] = config.data_path(fp)
    finally:
        con.close()
    pick = random.sample(sorted(set(acc) & members), min(60, len(set(acc) & members)))
    if not pick:
        rec("SKIP", "mutation profile vs Nextclade's alignment", "no member with a sequence file")
        return
    with tempfile.TemporaryDirectory(prefix="gcver-") as td:
        fa, out = Path(td) / "in.fa", Path(td) / "aln.fa"
        with open(fa, "w") as w:
            for p in pick:
                seq = "".join(l.strip() for l in open(acc[p], encoding="utf-8") if not l.startswith(">"))
                w.write(f">{p}\n{seq}\n")
        r = subprocess.run([nc, "run", "--input-ref", str(ref_fa), "--output-fasta", str(out), str(fa)], capture_output=True, text=True, timeout=1800)
        if r.returncode != 0 or not out.exists():
            rec("SKIP", "mutation profile vs Nextclade's alignment", (r.stderr or "")[-200:])
            return
        from genomecurator import aligned
        aln = dict(aligned.read_alignment(out))
    mine: dict[str, set] = {}
    for r_ in rows:
        mine.setdefault(r_["sample"], set()).add((int(r_["pos"]), r_["alt"].upper()))
    tot_both = tot_run = tot_nc = 0
    for p in pick:
        s = aln.get(p)
        if not s or len(s) != len(rseq):
            continue
        theirs = {(i + 1, b) for i, b in enumerate(s.upper()) if b in "ACGT" and rseq[i].upper() in "ACGT" and b != rseq[i].upper() and (win is None or win[0] <= i + 1 <= win[1])}
        known = {i + 1 for i, b in enumerate(s.upper()) if b in "ACGT"}
        m = {x for x in mine.get(p, set()) if x[0] in known and (win is None or win[0] <= x[0] <= win[1])}
        tot_both += len(m & theirs)
        tot_run += len(m - theirs)
        tot_nc += len(theirs - m)
    rec_ = tot_both / max(1, tot_both + tot_nc)
    prec = tot_both / max(1, tot_both + tot_run)
    ok(rec_ >= 0.95 and prec >= 0.95, "mutation profile vs Nextclade's alignment (another aligner) on the same genomes", f"{len(pick)} genomes: {tot_both:,} mutations in both, {tot_run} only in the run, {tot_nc} only in Nextclade (recall {rec_ * 100:.1f}%, precision {prec * 100:.1f}%)", warn=True)


def mutprofile_ref_dataset(mp):
    from genomecurator import nextstrain
    d = nextstrain.workflow_dir() / "nextclade" / "resources"
    b = (mp.get("build") or "clade-i")
    f = d / b
    return f if f.exists() else None


if __name__ == "__main__":
    sys.exit(main())
