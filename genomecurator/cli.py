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
