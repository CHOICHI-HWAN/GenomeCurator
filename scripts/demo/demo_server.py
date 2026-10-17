"""Demo-only server: the normal app, except that the schema check sees a *pretend* source change while the flag file exists.
Used by `make_gifs.py schema` to show what "update needed" looks like. Never used by the real app."""
import copy
import os
import sys
from pathlib import Path

import uvicorn

from genomecurator import schema_check as sc
from genomecurator.web import create_app

FLAG = Path(os.environ["GC_DEMO_FLAG"])
_n, _p = sc._ncbi_sample, sc._ppx_sample


def ncbi(c, taxon, n=sc.SAMPLE):
    rows, shape = _n(c, taxon, n)
    if FLAG.exists():   # pretend: NCBI dropped `host` and renamed `length` to `seqLength`
        rows = copy.deepcopy(rows)
        for r in rows:
            r.pop("host", None)
            r["seqLength"] = r.pop("length", None)
    return rows, shape


def ppx(c, slug, n=sc.SAMPLE):
    rows, shape = _p(c, slug, n)
    return (rows, {"top": ["results", "info"], "status": 200}) if FLAG.exists() else (rows, shape)   # pretend: data -> results


sc._ncbi_sample, sc._ppx_sample = ncbi, ppx
uvicorn.run(create_app(scheduler=False), host="127.0.0.1", port=int(sys.argv[1]), log_level="warning")
