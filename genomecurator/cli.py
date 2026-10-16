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


def main(argv=None):
    _utf8_console()
    p = argparse.ArgumentParser(prog="genomecurator")
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
    s = sub.add_parser("serve", help="run the web UI")
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--port", type=int, default=8000)
    s.add_argument("--no-scheduler", action="store_true")
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
    elif args.cmd == "dry-run":
        print(json.dumps(curate.dryrun(args.changeset), indent=2, ensure_ascii=False))
    elif args.cmd == "apply":
        if args.approve_all:
            ids = [r[0] for r in connect().execute(
                "SELECT id FROM changeset_item WHERE changeset_id=? AND decision='pending'", (args.changeset,))]
            curate.decide(ids, "approved")
        print(curate.apply(args.changeset))
    elif args.cmd == "serve":
        import uvicorn
        from .web import create_app
        if args.open:
            threading.Timer(1.5, lambda: webbrowser.open(f"http://{'127.0.0.1' if args.host in ('0.0.0.0', '::') else args.host}:{args.port}/")).start()
        uvicorn.run(create_app(scheduler=not args.no_scheduler), host=args.host, port=args.port)


if __name__ == "__main__":
    main()
