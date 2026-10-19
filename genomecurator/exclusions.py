"""Leave genomes out of the final analysis dataset - without deleting anything.

An exclusion is only a mark (who, why, when) kept in `dataset_exclusion`; the genome row, its sequence, its history and its QC all stay, and the
mark can be lifted at any time. Exports (`export.py`) skip marked genomes and record them in `excluded.tsv` / the manifest, so a dataset always
says what was left out and why. It is separate from `Removed` (a genome that NCBI dropped or that fell below the coverage floor): that is
about the curated set itself, this is a per-dataset decision. Works for public genomes (by accession) and for the private store (by id)."""
from datetime import datetime

from .db import audit, connect

SCOPES = ("public", "private")


def _exists(con, scope: str, watch_id: int, i: str) -> bool:
    if scope == "private":
        return con.execute("SELECT 1 FROM private_genome WHERE watch_id=? AND id=?", (watch_id, i)).fetchone() is not None
    return con.execute("SELECT 1 FROM genome WHERE watch_id=? AND accession=?", (watch_id, i)).fetchone() is not None


def add(watch_id: int, ids: list[str], scope: str = "public", reason: str = "") -> dict:
    """Mark genomes as excluded (idempotent: marking again only updates the reason). Unknown ids are reported, not stored."""
    if scope not in SCOPES:
        raise ValueError("scope는 public 또는 private입니다.")
    con = connect()
    try:
        if not con.execute("SELECT 1 FROM watchlist WHERE id=?", (watch_id,)).fetchone():
            raise ValueError("바이러스를 찾을 수 없습니다.")
        now = datetime.now().isoformat(timespec="seconds")
        done, new, unknown = 0, 0, []
        for i in dict.fromkeys(str(x).strip() for x in ids if str(x).strip()):
            if not _exists(con, scope, watch_id, i):
                unknown.append(i)
                continue
            old = con.execute("SELECT 1 FROM dataset_exclusion WHERE scope=? AND watch_id=? AND id=?", (scope, watch_id, i)).fetchone()
            con.execute("INSERT OR REPLACE INTO dataset_exclusion(scope,watch_id,id,reason,excluded_at) VALUES(?,?,?,?,?)", (scope, watch_id, i, reason.strip(), now))
            done += 1
            new += 0 if old else 1
        if done:
            audit(con, "dataset_exclude", f"{scope} 바이러스 {watch_id}: {new}건 제외 표시" + (f" ({reason.strip()})" if reason.strip() else ""))
        con.commit()
        return {"excluded": done, "new": new, "unknown": unknown}
    finally:
        con.close()


def remove(watch_id: int, ids: list[str], scope: str = "public") -> dict:
    """Lift the mark (the genomes were never deleted)."""
    con = connect()
    try:
        n = 0
        for i in dict.fromkeys(str(x).strip() for x in ids if str(x).strip()):
            n += con.execute("DELETE FROM dataset_exclusion WHERE scope=? AND watch_id=? AND id=?", (scope, watch_id, i)).rowcount
        if n:
            audit(con, "dataset_include", f"{scope} 바이러스 {watch_id}: {n}건 제외 해제")
        con.commit()
        return {"included": n}
    finally:
        con.close()


def listing(watch_id: int, scope: str | None = None) -> list[dict]:
    con = connect()
    try:
        q, a = "SELECT scope, id, reason, excluded_at FROM dataset_exclusion WHERE watch_id=?", [watch_id]
        if scope:
            q += " AND scope=?"
            a.append(scope)
        return [dict(r) for r in con.execute(q + " ORDER BY scope, id", a)]
    finally:
        con.close()


def marked(watch_id: int, scope: str = "public") -> dict[str, dict]:
    """{id: {reason, excluded_at}} for quick lookups."""
    return {r["id"]: r for r in listing(watch_id, scope)}
