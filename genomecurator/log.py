"""Logging for GenomeCurator: one rotating file (data/logs/genomecurator.log) plus the console.

Every line carries a context (job / watch / changeset) so a problem can be followed across modules:
    2026-10-17 16:02:11,402 INFO  genomecurator.sources [job=ab12cd34 watch=Monkeypox virus] GET ... -> 200 1.2s
Use `get("name")` for a logger and `with context(job=...)` to tag everything done inside it. Levels: DEBUG every HTTP call,
INFO milestones (start/finish/counts), WARNING a retry or a degraded fallback, ERROR a failed step (with traceback).
GC_LOG_LEVEL changes the file level (default INFO); GC_LOG_CONSOLE_LEVEL the console (default WARNING)."""
import contextvars
import logging
import logging.handlers
import os
import re
import sys
from contextlib import contextmanager

from . import config

_ctx: contextvars.ContextVar[dict] = contextvars.ContextVar("gc_log_ctx", default={})
ROOT = "genomecurator"
logging.getLogger(ROOT).addHandler(logging.NullHandler())
_ready = False


class _Context(logging.Filter):
    def filter(self, record):
        c = _ctx.get()
        record.ctx = (" [" + " ".join(f"{k}={v}" for k, v in c.items()) + "]") if c else ""
        return True


def get(name: str) -> logging.Logger:
    return logging.getLogger(f"{ROOT}.{name}")


@contextmanager
def context(**kv):
    """Tag every log line written inside the block (and by code it calls) with key=value pairs."""
    token = _ctx.set({**_ctx.get(), **{k: v for k, v in kv.items() if v not in (None, "")}})
    try:
        yield
    finally:
        _ctx.reset(token)


def set_context(**kv):
    """For threads that start without the parent's context (jobs): set it once at the top of the thread."""
    _ctx.set({k: v for k, v in kv.items() if v not in (None, "")})


def log_dir():
    return config.DATA_DIR / "logs"


def setup(console: bool = True) -> None:
    """Idempotent. Called by the CLI and the web app; tests leave logging unconfigured."""
    global _ready
    if _ready:
        return
    root = logging.getLogger(ROOT)
    root.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s %(levelname)-5s %(name)s%(ctx)s %(message)s")
    flt = _Context()
    try:
        log_dir().mkdir(parents=True, exist_ok=True)
        fh = logging.handlers.RotatingFileHandler(log_dir() / "genomecurator.log", maxBytes=5_000_000, backupCount=5, encoding="utf-8")
        fh.setLevel(os.environ.get("GC_LOG_LEVEL", "INFO").upper())
        fh.setFormatter(fmt)
        fh.addFilter(flt)
        root.addHandler(fh)
    except OSError as e:   # a read-only data folder must not stop the app
        print(f"warning: cannot write log file: {e}", file=sys.stderr)
    if console:
        ch = logging.StreamHandler(sys.stderr)
        ch.setLevel(os.environ.get("GC_LOG_CONSOLE_LEVEL", "WARNING").upper())
        ch.setFormatter(fmt)
        ch.addFilter(flt)
        root.addHandler(ch)
    _ready = True
    from . import __version__
    get("log").info("logging started: Virus Genome Curator %s (file level %s) -> %s", __version__, os.environ.get("GC_LOG_LEVEL", "INFO").upper(), log_dir())


def tail(n: int = 200, level: str = "", contains: str = "") -> list[str]:
    """Last n lines of the current log file (optionally only WARNING/ERROR), for the web log viewer."""
    try:
        lines = (log_dir() / "genomecurator.log").read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    if contains:      # e.g. "job=ab12cd34": every line of one job, plus the traceback lines that follow a matching line
        out, keep = [], False
        for ln in lines:
            starts = bool(re.match(r"\d{4}-\d{2}-\d{2} ", ln))
            keep = (contains in ln) if starts else keep
            if keep:
                out.append(ln)
        lines = out
    if level:
        keep = {"WARNING": ("WARNING", "ERROR", "CRITICAL"), "ERROR": ("ERROR", "CRITICAL")}.get(level.upper())
        if keep:
            lines = [ln for ln in lines if any(f" {k} " in ln[:40] or f" {k}" in ln[:30] for k in keep) or ln.startswith(("Traceback", "  "))]
    return lines[-n:]
