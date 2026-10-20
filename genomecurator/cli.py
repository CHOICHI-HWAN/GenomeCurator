import argparse
import json
import sys
import threading
import webbrowser

from . import config, curate
from .db import connect


def _utf8_console():
    """Korean text must not crash print() on consoles with a legacy code page (e.g. Windows cp949)."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass


def _main(argv=None):
    _utf8_console()
    from . import log
    log.setup()
    from . import DB_SCHEMA_VERSION, __version__
    p = argparse.ArgumentParser(prog="genomecurator")
    p.add_argument("--version", action="version", version=f"genomecurator {__version__} (Virus Genome Metadata Curator, DB schema {DB_SCHEMA_VERSION})")
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("watch-add", help="register a virus (taxon name or id)")
    a.add_argument("name")
    a.add_argument("--taxon", help="NCBI taxon name/id (default: name)")
    a.add_argument("--filter", action="append", default=[], metavar="K=V",
                   help="NCBI filter, e.g. filter.complete_only=true")
    sub.add_parser("watch-list")
    c = sub.add_parser("check", help="fetch sources and create a changeset (no local data change)")
    c.add_argument("name")
    d = sub.add_parser("dry-run", help="validate a changeset's download plan")
    d.add_argument("changeset", type=int)
    ap = sub.add_parser("apply", help="apply approved items of a READY changeset")
    ap.add_argument("changeset", type=int)
    ap.add_argument("--approve-all", action="store_true", help="CLI shortcut; the web page is the intended review path")
    sc = sub.add_parser("schema-check", help="dry-run: compare the live schema of NCBI/BV-BRC/Pathoplexus and the local DB with what the code expects")
    sc.add_argument("--daily", action="store_true", help="the light per-virus check (3 records per source) that the scheduler runs every day")
    rc = sub.add_parser("reconcile", help="show the cross-source agreement report stored by the last check")
    rc.add_argument("name")
    vf = sub.add_parser("verify", help="integrity checks on local data (FASTA length/md5, dates, orphans)")
    vf.add_argument("--no-files", action="store_true", help="skip reading FASTA files")
    bk = sub.add_parser("backup", help="write a portable backup (curated data + FASTA) to a .zip")
    bk.add_argument("path", nargs="?", default="genomecurator-backup.zip")
    bk.add_argument("--with-private", action="store_true", help="also include the private store; such a backup must not be shared")
    ex = sub.add_parser("export", help="write metadata.tsv + sequences.fasta (+ manifest) for other programs, e.g. the Nextstrain mpox workflow")
    ex.add_argument("name", help="the registered virus")
    ex.add_argument("--profile", default="nextstrain", help="nextstrain (default) or plain")
    ex.add_argument("--out", help="folder to write (default: <data>/exports/<virus>/<timestamp>)")
    ex.add_argument("--nextstrain-dir", help="the phylogenetic/ folder of a nextstrain/mpox checkout: write the files into <it>/data/genomecurator/ and a config overlay into <it>/build-configs/genomecurator/")
    ex.add_argument("--min-length", type=int)
    ex.add_argument("--completeness", help="comma list, e.g. 'Complete,Near complete'")
    ex.add_argument("--exclude-qc", help="comma list of QC statuses to leave out, e.g. 'bad,error'")
    ex.add_argument("--clade")
    ex.add_argument("--with-private", action="store_true", help="also write the private store as metadata-private.tsv / sequences-private.fasta (never mixed into the public files)")
    ex.add_argument("--with-missing-sequences", action="store_true", help="also list genomes that have no sequence file (metadata only)")
    pc = sub.add_parser("private-check", help="dry run: validate a private metadata TSV + FASTA (nothing is stored)")
    pi = sub.add_parser("private-import", help="import private genomes (metadata TSV + FASTA) into the separate private store")
    for x in (pc, pi):
        x.add_argument("name", help="the registered virus these genomes belong to")
        x.add_argument("--metadata", required=True, help="TSV with an id column (PPX_accession or accession), strain, date, country, ...")
        x.add_argument("--fasta", required=True, help="FASTA whose ids are those ids")
    dc = sub.add_parser("dataset-create", help="save a named dataset (rules over the curated genomes) - e.g. only clade Ib")
    dc.add_argument("name", help="the registered virus")
    dc.add_argument("dataset", help="the dataset's name")
    for opt in ("clade", "lineage", "outbreak", "region", "country", "host", "completeness"):
        dc.add_argument(f"--{opt}", help="comma list, e.g. --clade Ib or --clade Ia,Ib")
    dc.add_argument("--exclude-qc", help="comma list of QC statuses to leave out, e.g. bad,error")
    dc.add_argument("--year-min", type=int)
    dc.add_argument("--year-max", type=int)
    dc.add_argument("--min-length", type=int)
    dc.add_argument("--max-length", type=int)
    dc.add_argument("--include", help="comma list of accessions that always belong")
    dc.add_argument("--leave-out", help="comma list of accessions that never belong")
    dc.add_argument("--private", action="store_true", help="the dataset also contains the private store")
    dc.add_argument("--with-missing-sequences", action="store_true", help="also members without a sequence file")
    dc.add_argument("--build", choices=("mpxv", "hmpxv1", "clade-i"), help="the nextstrain/mpox build the exported config suggests")
    dc.add_argument("--description", default="")
    dc.add_argument("--update", action="store_true", help="replace an existing dataset of that name")
    dl = sub.add_parser("dataset-list", help="list the saved datasets of a virus")
    dl.add_argument("name")
    de = sub.add_parser("dataset-export", help="write a saved dataset (metadata.tsv, sequences.fasta, manifest ...)")
    de.add_argument("name")
    de.add_argument("dataset")
    de.add_argument("--out")
    de.add_argument("--nextstrain-dir")
    dx = sub.add_parser("dataset-delete", help="delete a dataset's definition (genomes and exported files stay)")
    dx.add_argument("name")
    dx.add_argument("dataset")
    na = sub.add_parser("analyze", help="run the Nextstrain mpox workflow on a saved dataset (every run is kept); `analyze <virus> --status` shows the setup")
    na.add_argument("name", help="the registered virus")
    na.add_argument("dataset", nargs="?", help="the dataset's name")
    na.add_argument("--status", action="store_true", help="environment, workflow checkout and viewer")
    na.add_argument("--clone", action="store_true", help="get / update the nextstrain/mpox workflow checkout")
    na.add_argument("--list", action="store_true", help="the kept runs of the dataset")
    na.add_argument("--runtime", choices=("ambient", "docker"), default="ambient")
    na.add_argument("--cores", type=int)
    na.add_argument("--force", action="store_true", help="run again even when members, settings and workflow are the same as a finished run")
    na.add_argument("--note", default="")
    pp = sub.add_parser("pathoplexus", help="prepare private genomes for submission to Pathoplexus (checks + files; nothing is uploaded)")
    pp.add_argument("name", help="the registered virus")
    pp.add_argument("--check", action="store_true", help="show what is wrong with each private genome and stop")
    pp.add_argument("--ids", help="comma list of private ids (default: all that have no errors)")
    pp.add_argument("--terms", default="OPEN", choices=("OPEN", "RESTRICTED"))
    pp.add_argument("--until", default="", help="YYYY-MM-DD, for RESTRICTED (at most one year from today)")
    pp.add_argument("--group-id", default="", help="your Pathoplexus group id (written into submit.sh)")
    pp.add_argument("--list", action="store_true", help="the packages prepared so far")
    vc = sub.add_parser("vocab", help="what the term dictionary would change in a virus's curated data (nothing is written)")
    vc.add_argument("name", help="the registered virus")
    vc.add_argument("--all", action="store_true", help="also list values that are already standard")
    mk = sub.add_parser("mask-consensus", help="columns where Squirrel's and Nextclade's alignments disagree -> a BED for the Nextstrain workflow's mask.maskfile (nothing else is changed)")
    mk.add_argument("squirrel_aln")
    mk.add_argument("squirrel_ref", help="Squirrel's reference FASTA (PX667572 / NC_063383)")
    mk.add_argument("nextclade_aln", help="Nextclade --output-fasta of the same genomes")
    mk.add_argument("target_ref", help="the reference Nextclade aligned to (reference.fasta of its dataset)")
    mk.add_argument("--min-genomes", type=int, default=3, help="mask a column when at least this many genomes disagree there (default 3)")
    mk.add_argument("-o", "--out", default="", help="BED file to write (default: only a summary)")
    xc = sub.add_parser("exclude", help="leave genomes out of the final analysis dataset (nothing is deleted); --undo puts them back")
    xc.add_argument("name", help="the registered virus")
    xc.add_argument("ids", nargs="*", help="accessions (or private ids with --private)")
    xc.add_argument("--reason", default="", help="why (recorded in excluded.tsv)")
    xc.add_argument("--private", action="store_true", help="the ids are private-store ids")
    xc.add_argument("--undo", action="store_true", help="lift the exclusion")
    xc.add_argument("--list", action="store_true", help="show what is excluded")
    pl = sub.add_parser("private-list", help="list the private genomes of a virus")
    pl.add_argument("name")
    q = sub.add_parser("qc", help="sequence statistics and Nextclade for one virus")
    q.add_argument("name")
    q.add_argument("--private", action="store_true", help="run on the private store instead of the public genomes")
    q.add_argument("--force", action="store_true", help="recompute even when nothing changed")
    rs = sub.add_parser("restore", help="restore a backup .zip into the data folder")
    rs.add_argument("path")
    rs.add_argument("--force", action="store_true", help="replace an existing database")
    s = sub.add_parser("serve", help="run the web UI")
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--port", type=int, default=8000)
    s.add_argument("--no-scheduler", action="store_true")
    s.add_argument("--superset", action="store_true", help="also start Apache Superset (port 8089) next to the app")
    s.add_argument("--open", action="store_true", help="open the browser when the server is up")
    args = p.parse_args(argv)

    if args.cmd == "watch-add":
        filters = dict(f.split("=", 1) for f in args.filter)
        print("watch id", curate.add_watch(args.name, args.taxon, filters))
    elif args.cmd == "watch-list":
        for r in connect().execute("SELECT id,name,taxon,taxon_id,last_checked FROM watchlist"):
            print(tuple(r))
    elif args.cmd == "check":
        cs = curate.check(args.name)
        print(f"changeset {cs}" if cs else "no changes")
        if cs:
            print(json.dumps(curate.dryrun(cs), indent=2, ensure_ascii=False))
    elif args.cmd == "schema-check":
        from . import schema_check
        r = schema_check.run_light() if args.daily else schema_check.run()
        word = {"ok": "소스 코드 업데이트 필요 없음", "review": "확인 권장", "update": "소스 코드 업데이트 필요"}[r["verdict"]]
        print(f"[{r['verdict'].upper()}] {word} (fail {r['fail']}, warn {r['warn']}, info {r['info']}) - {('per-virus daily check' if r.get('mode') == 'daily' else 'sample from ' + r['watch'])}")
        for f in r["findings"]:
            if f["severity"] in ("fail", "warn"):
                print(f"  {f['severity'].upper():4} {f.get('virus', ''):24} {f['source']:11} {f['kind']:10} {f['path']}: {f['detail']} -> {f['hint']}")
        raise SystemExit(1 if r["verdict"] == "update" else 0)
    elif args.cmd == "reconcile":
        from . import reconcile
        w = connect().execute("SELECT id FROM watchlist WHERE name=?", (args.name,)).fetchone()
        r = reconcile.load(w[0]) if w else None
        if not r:
            raise SystemExit("no report yet: run check first")
        print(json.dumps(r, indent=2, ensure_ascii=False))
    elif args.cmd == "verify":
        from . import verify
        bad = 0
        for c in verify.run(check_files=not args.no_files):
            mark = {"ok": "OK  ", "warn": "WARN", "fail": "FAIL"}[c["status"]]
            print(f"[{mark}] {c['name']}: {c['count']} {c['detail']} {c['examples'] or ''}")
            bad += c["status"] == "fail"
        raise SystemExit(1 if bad else 0)
    elif args.cmd == "dry-run":
        print(json.dumps(curate.dryrun(args.changeset), indent=2, ensure_ascii=False))
    elif args.cmd == "apply":
        if args.approve_all:
            ids = [r[0] for r in connect().execute(
                "SELECT id FROM changeset_item WHERE changeset_id=? AND decision='pending'", (args.changeset,))]
            curate.decide(ids, "approved")
        print(curate.apply(args.changeset))
    elif args.cmd in ("dataset-create", "dataset-list", "dataset-export", "dataset-delete"):
        from pathlib import Path
        from . import datasets, export
        con = connect()
        try:
            w = con.execute("SELECT id FROM watchlist WHERE name=?", (args.name,)).fetchone()
        finally:
            con.close()
        if not w:
            raise SystemExit(f"unknown watch: {args.name}")
        try:
            if args.cmd == "dataset-list":
                for d in datasets.listing(w["id"]):
                    print(f"{d['name']}\t{d['summary']}\t{'비공개 포함' if d['include_private'] else ''}\t{d['last_export'] or ''}")
            elif args.cmd == "dataset-create":
                csv_ = lambda s: [x.strip() for x in s.split(",") if x.strip()] if s else []
                rules = {k: csv_(getattr(args, k)) for k in ("clade", "lineage", "outbreak", "region", "country", "host", "completeness")}
                rules.update(qc_not=csv_(args.exclude_qc), year_min=args.year_min, year_max=args.year_max, min_length=args.min_length, max_length=args.max_length,
                             require_sequence=not args.with_missing_sequences)
                old = datasets.find(w["id"], args.dataset) if args.update else None
                d = datasets.save(w["id"], args.dataset, rules, args.description, args.private, csv_(args.include), csv_(args.leave_out), args.build or "", old["id"] if old else None)
                print(json.dumps({"id": d["id"], "name": d["name"], "rules": datasets.describe(d["rules"])}, ensure_ascii=False))
            else:
                d = datasets.find(w["id"], args.dataset)
                if not d:
                    raise ValueError(f"데이터셋이 없습니다: {args.dataset}")
                if args.cmd == "dataset-delete":
                    print(json.dumps({"deleted": datasets.delete(d["id"])}))
                else:
                    pr = lambda pct, msg: print(f"[{round(pct):3d}%] {msg}", file=sys.stderr)
                    print(json.dumps(export.write(w["id"], Path(args.out) if args.out else None, "nextstrain", pr, nextstrain_dir=args.nextstrain_dir, dataset=d), ensure_ascii=False, indent=2))
        except ValueError as e:
            raise SystemExit(str(e))
    elif args.cmd == "analyze":
        from . import datasets, nextstrain
        pr = lambda pct, msg: print(f"[{round(pct):3d}%] {msg}", file=sys.stderr)
        con = connect()
        try:
            w = con.execute("SELECT id FROM watchlist WHERE name=?", (args.name,)).fetchone()
        finally:
            con.close()
        if not w:
            raise SystemExit(f"unknown watch: {args.name}")
        try:
            if args.status:
                print(json.dumps(nextstrain.status(), ensure_ascii=False, indent=2))
            elif args.clone:
                print(json.dumps(nextstrain.clone_workflow(pr), ensure_ascii=False))
            else:
                d = datasets.find(w["id"], args.dataset or "")
                if not d:
                    raise ValueError(f"데이터셋이 없습니다: {args.dataset}")
                if args.list:
                    for r in nextstrain.list_runs(d["id"]):
                        print("\t".join(str(r[k] or "") for k in ("name", "status", "started", "n_public", "n_private", "workflow_commit", "note")))
                else:
                    r = nextstrain.run(d, args.runtime, args.cores, args.force, args.note, pr)
                    print(json.dumps({k: r[k] for k in ("id", "name", "status", "auspice_file", "run_dir") if k in r} | ({"reused": True} if r.get("reused") else {}), ensure_ascii=False, indent=2))
        except (ValueError, RuntimeError) as e:
            raise SystemExit(str(e))
    elif args.cmd == "pathoplexus":
        from . import pathoplexus
        con = connect()
        try:
            w = con.execute("SELECT id FROM watchlist WHERE name=?", (args.name,)).fetchone()
        finally:
            con.close()
        if not w:
            raise SystemExit(f"unknown watch: {args.name}")
        try:
            if args.list:
                for p in pathoplexus.list_packages(w["id"]):
                    print("\t".join(str(p[k] or "") for k in ("id", "status", "created", "n", "terms", "dir")))
            else:
                cand = pathoplexus.candidates(w["id"])
                if args.check or not cand["rows"]:
                    print(f"organism: {cand['organism'] or '(unknown)'}  required: {', '.join(cand['required'])}")
                    for r in cand["rows"]:
                        print(f"{r['status'].upper():5} {r['id']}" + "".join(f"\n        ! {m}" for m in r["errors"]) + "".join(f"\n        ~ {m}" for m in r["warnings"]))
                    if not cand["rows"]:
                        print("비공개 저장소가 비어 있습니다.")
                else:
                    ids = [i for i in (args.ids or "").split(",") if i] or [r["id"] for r in cand["rows"] if not r["errors"]]
                    print(json.dumps(pathoplexus.build_package(w["id"], ids, args.terms, args.until, args.group_id), ensure_ascii=False, indent=2))
        except ValueError as e:
            raise SystemExit(str(e))
    elif args.cmd == "mask-consensus":
        from . import aligned, mutprofile
        dis = aligned.disagreement(args.squirrel_aln, args.squirrel_ref, args.nextclade_aln, args.target_ref)
        runs = aligned.consensus_mask(dis, args.min_genomes)
        n = sum(b - a + 1 for a, b in runs)
        worst = sorted(dis["per_genome"].items(), key=lambda x: -x[1])[:5]
        print(f"{dis['n']}개 유전체 비교, 참조 길이 {dis['length']:,} — 마스크 {len(runs)}구간 {n}열 ({100 * n / dis['length']:.3f}%), 불일치가 가장 많은 유전체: " + ", ".join(f"{k}({v})" for k, v in worst))
        if args.out:
            seqid = mutprofile.read_fasta(__import__("pathlib").Path(args.target_ref))[0].split()[0]
            __import__("pathlib").Path(args.out).write_text(aligned.mask_bed(runs, seqid), encoding="utf-8")
            print("BED 저장:", args.out)
    elif args.cmd == "vocab":
        from . import vocab
        con = connect()
        try:
            w = con.execute("SELECT id, name, taxon FROM watchlist WHERE name=?", (args.name,)).fetchone()
            if not w:
                raise SystemExit(f"unknown watch: {args.name}")
            pv = vocab.preview(con, w["id"], vocab.scope_for(w["name"], w["taxon"] or ""))
        finally:
            con.close()
        for f, d in pv.items():
            print(f"{f}: {d['genomes']:,} genomes, {d['distinct']} distinct -> {d['distinct_after']} | 바뀜 {d['changes']:,} | 검토 필요 {d['review']:,}")
            for it in d["items"]:
                if args.all or it["status"] not in ("same", "nodict"):
                    print(f"   {it['status']:<10} {it['raw']!r} -> {it['standard']!r}  ({it['n']:,})  {it['rule']}")
    elif args.cmd == "exclude":
        from . import exclusions
        con = connect()
        try:
            w = con.execute("SELECT id FROM watchlist WHERE name=?", (args.name,)).fetchone()
        finally:
            con.close()
        if not w:
            raise SystemExit(f"unknown watch: {args.name}")
        scope = "private" if args.private else "public"
        if args.list or not args.ids:
            for e in exclusions.listing(w["id"], None if args.list and not args.private else scope):
                print("\t".join(str(e[k] or "") for k in ("scope", "id", "excluded_at", "reason")))
        elif args.undo:
            print(json.dumps(exclusions.remove(w["id"], args.ids, scope), ensure_ascii=False))
        else:
            print(json.dumps(exclusions.add(w["id"], args.ids, scope, args.reason), ensure_ascii=False))
    elif args.cmd in ("private-check", "private-import", "private-list"):
        from . import private
        from pathlib import Path
        con = connect()
        try:
            w = con.execute("SELECT id FROM watchlist WHERE name=?", (args.name,)).fetchone()
        finally:
            con.close()
        if not w:
            raise SystemExit(f"unknown watch: {args.name}")
        if args.cmd == "private-list":
            for r in private.list_rows(w["id"]):
                print("\t".join(str(r[k] or "") for k in ("id", "strain", "date", "country", "clade", "length", "qc_stats", "qc_nextclade")))
        else:
            mt, fa = Path(args.metadata).read_text(encoding="utf-8"), Path(args.fasta).read_text(encoding="utf-8")
            if args.cmd == "private-check":
                rep = private.check(w["id"], mt, fa)
                print(json.dumps({k: rep[k] for k in ("ok", "n", "new", "update", "id_column", "errors", "warnings")}, ensure_ascii=False, indent=2))
                if not rep["ok"]:
                    raise SystemExit(1)
            else:
                try:
                    print(json.dumps(private.import_texts(w["id"], mt, fa, Path(args.metadata).name), ensure_ascii=False, indent=2))
                except ValueError as e:
                    raise SystemExit(str(e))
    elif args.cmd in ("export", "qc"):
        from pathlib import Path
        con = connect()
        try:
            w = con.execute("SELECT id FROM watchlist WHERE name=?", (args.name,)).fetchone()
        finally:
            con.close()
        if not w:
            raise SystemExit(f"unknown watch: {args.name}")
        pr = lambda pct, msg: print(f"[{round(pct):3d}%] {msg}", file=sys.stderr)
        if args.cmd == "qc":
            from . import qc
            print(json.dumps((qc.run_private if args.private else qc.run)(w["id"], pr, args.force), ensure_ascii=False, indent=2))
        else:
            from . import export
            csv_ = lambda s: [x.strip() for x in s.split(",") if x.strip()] if s else None
            r = export.write(w["id"], Path(args.out) if args.out else None, args.profile, pr, nextstrain_dir=args.nextstrain_dir, include_private=args.with_private, require_sequence=not args.with_missing_sequences,
                             min_length=args.min_length, completeness_in=csv_(args.completeness), qc_not=csv_(args.exclude_qc), clade=args.clade)
            print(json.dumps(r, ensure_ascii=False, indent=2))
    elif args.cmd == "backup":
        from . import backup
        r = backup.export_backup(args.path, with_private=args.with_private)
        print(f"backup written: {r['path']} ({r['size_mb']} MB) - tables {r['tables']}, FASTA files {r['fasta_files']}")
    elif args.cmd == "restore":
        from . import backup
        try:
            r = backup.restore_backup(args.path, force=args.force)
        except (FileExistsError, ValueError, RuntimeError) as e:
            raise SystemExit(str(e))
        print(f"restored backup from {r['created']}: tables {r['tables']}, FASTA files {r['fasta_files']} -> {config.DATA_DIR}")
    elif args.cmd == "serve":
        import uvicorn
        from .web import create_app
        if args.open:
            threading.Timer(1.5, lambda: webbrowser.open(f"http://{'127.0.0.1' if args.host in ('0.0.0.0', '::') else args.host}:{args.port}/")).start()
        if args.superset:
            from . import analysis
            threading.Thread(target=lambda: _start_superset(analysis), daemon=True).start()
        uvicorn.run(create_app(scheduler=not args.no_scheduler), host=args.host, port=args.port)


def _start_superset(analysis):
    try:
        print("Superset:", analysis.start())
    except Exception as e:
        print("Superset을 시작하지 못했습니다:", e)


def main(argv=None):
    from . import log, sources
    from .db import audit_now
    try:
        return _main(argv)
    except sources.SourceUnavailable as e:   # a source failed: nothing was changed; say so plainly instead of a traceback
        log.get("cli").error("command failed: %s", e)
        audit_now("command_failed", str(e))
        raise SystemExit(f"오류: {e}\n(자세한 기록: {log.log_dir() / 'genomecurator.log'})") from e


if __name__ == "__main__":
    main()
