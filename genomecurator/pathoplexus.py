"""Prepare private genomes for submission to Pathoplexus - and nothing more: the app never uploads and never sees a password or token.

Submitting publishes data (open data is released to INSDC, and restricted data can only be shortened or opened later, never extended), so the last step is always
the person's own: the files this module writes are uploaded on pathoplexus.org (which has a review and an approval step before anything is released) or with the
generated `submit.sh`, which reads the token from an environment variable of the person's own shell.

What is prepared (following Pathoplexus' documentation: a metadata TSV whose `id` column equals the FASTA ids, with the organism's own metadata template):
* the organism's *template* is read from the site (`/<organism>/submission/template`, a public GET) so the columns follow whatever Pathoplexus asks for today;
* private metadata is mapped to it (country, collection date as YYYY / YYYY-MM / YYYY-MM-DD, authors as `Last, First; ...`, affiliations, host, admin areas, sample id; any extra
  column of the private sheet that is named like a template column - sequencingInstrument, sequencingProtocol, depthOfCoverage, comment ... - is passed through);
* every genome is checked before anything is written: required fields, formats, authors, host, whether the same sequence is already public, QC status;
* a package folder with `metadata.tsv`, `sequences.fasta`, `submission.json` (what, which versions, SHA-256), `submit.sh` and a README with the steps;
* afterwards the person can record that it was submitted and the accessions Pathoplexus returned."""
import hashlib
import json
import re
import urllib.request
from datetime import date, datetime, timedelta
from pathlib import Path

from . import config, export, log, private
from .db import audit, connect

_log = log.get("pathoplexus")
SITE = "https://pathoplexus.org"
BACKEND = "https://backend.pathoplexus.org"
# the columns of the mpox template as of this writing; the live template replaces it when it can be read
FALLBACK_TEMPLATE = ["id", "geoLocCountry", "sampleCollectionDate", "authors", "authorAffiliations", "host", "cultureId", "depthOfCoverage", "geoLocAdmin1", "geoLocAdmin2",
                     "sequencingInstrument", "sequencingProtocol", "specimenCollectorSampleId", "versionComment", "comment", "files.rawReads"]
REQUIRED = ["id", "geoLocCountry", "sampleCollectionDate", "authors", "authorAffiliations", "host"]   # the template lists these first; Pathoplexus' own check has the last word
AUTHOR = re.compile(r"^\s*[A-Za-z][A-Za-z'\- ]*,\s*[A-Za-z.' \-]*\s*$")
PASSTHROUGH_SKIP = {"id", "PPX_accession", "accession", "strain", "date", "country", "division", "location", "region", "length", "host", "institution", "authors"}


def organism_for(watch_id: int) -> str:
    con = connect()
    try:
        r = con.execute("SELECT pathoplexus_organism, taxon FROM watchlist WHERE id=?", (watch_id,)).fetchone()
    finally:
        con.close()
    if r and r["pathoplexus_organism"]:
        return r["pathoplexus_organism"]
    return "mpox" if r and "monkeypox" in (r["taxon"] or "").lower() else ""


def packages_dir() -> Path:
    return config.DATA_DIR / "pathoplexus"


def template(organism: str, refresh: bool = False) -> list[str]:
    """The organism's metadata template columns: the cached copy, else the live one (public GET), else the built-in mpox list."""
    f = packages_dir() / f"template-{organism}.tsv"
    if f.exists() and not refresh:
        return f.read_text(encoding="utf-8").splitlines()[0].split("\t")
    try:
        req = urllib.request.Request(f"{SITE}/{organism}/submission/template?format=tsv", headers={"User-Agent": "GenomeCurator"})
        with urllib.request.urlopen(req, timeout=20) as r:
            text = r.read().decode("utf-8")
        head = text.splitlines()[0].split("\t")
        if "id" in head and len(head) > 3:
            f.parent.mkdir(parents=True, exist_ok=True)
            f.write_text(head and "\t".join(head) + "\n", encoding="utf-8")
            return head
    except Exception:
        _log.warning("could not read the Pathoplexus template for %s", organism, exc_info=True)
    return f.read_text(encoding="utf-8").splitlines()[0].split("\t") if f.exists() else list(FALLBACK_TEMPLATE)


def pp_date(d: str) -> str:
    """Our '2025-03-XX' / '2025-XX-XX' -> Pathoplexus' 2025-03 / 2025 (no X's)."""
    d = (d or "").strip()
    m = re.fullmatch(r"(\d{4})(?:-(\d{2}|XX))?(?:-(\d{2}|XX))?", d)
    if not m:
        return ""
    y, mo, dd = m.groups()
    if not mo or mo == "XX":
        return y
    if not dd or dd == "XX":
        return f"{y}-{mo}"
    return f"{y}-{mo}-{dd}"


def _check_authors(a: str) -> str:
    a = (a or "").strip()
    if not a:
        return "authors가 비어 있습니다 (Pathoplexus 필수: 'Smith, Anna; Perez, Tom J.;' 형식)"
    parts = [p for p in a.split(";") if p.strip()]
    bad = [p.strip() for p in parts if not AUTHOR.match(p)]
    if bad:
        return f"authors 형식이 맞지 않습니다 ('성, 이름; 성, 이름' 형식, 성은 알파벳): {', '.join(bad[:3])}"
    return ""


def map_row(rec: dict, tmpl: list[str]) -> dict:
    """One private genome (rows_for_export record) as a Pathoplexus metadata row. Only columns that exist in the template are produced."""
    m = rec["meta"]
    row = {c: "" for c in tmpl}
    put = lambda k, v: row.__setitem__(k, v) if k in row else None
    put("id", rec["id"])
    put("geoLocCountry", m.get("country", ""))
    put("sampleCollectionDate", pp_date(m.get("date", "")))
    put("authors", (m.get("authors") or "").strip())
    put("authorAffiliations", (m.get("institution") or m.get("authorAffiliations") or "").strip())
    put("host", m.get("host", ""))
    put("geoLocAdmin1", m.get("division", ""))
    put("geoLocAdmin2", m.get("location", ""))
    put("specimenCollectorSampleId", m.get("strain", ""))
    for k, v in m.items():                       # columns of the private sheet named like a template column are passed through
        if k in row and k not in PASSTHROUGH_SKIP and not row[k] and str(v).strip():
            row[k] = str(v).strip()
    return row


def candidates(watch_id: int) -> dict:
    """Every private genome of the virus with what is wrong (error) or doubtful (warning) about submitting it, plus the template that was used."""
    org = organism_for(watch_id)
    tmpl = template(org) if org else list(FALLBACK_TEMPLATE)
    con = connect()
    try:
        pubmd5 = {r[0]: r[1] for r in con.execute("SELECT s.md5, g.accession FROM genome_sequence s JOIN genome g ON g.genome_id=s.genome_id WHERE s.md5 IS NOT NULL")}
        pkg = {}
        for r in con.execute("SELECT id, status, ids_json, created FROM submission_package WHERE watch_id=? ORDER BY id", (watch_id,)):
            for i in json.loads(r["ids_json"] or "[]"):
                pkg[i] = {"package": r["id"], "status": r["status"], "created": r["created"]}
        md5 = {r["id"]: r["seq_md5"] for r in con.execute("SELECT id, seq_md5 FROM private_genome WHERE watch_id=?", (watch_id,))}
    finally:
        con.close()
    from . import exclusions
    ex = exclusions.marked(watch_id, "private")
    out = []
    for rec in private.rows_for_export(watch_id):
        row = map_row(rec, tmpl)
        errors, warns = [], []
        for c in REQUIRED:
            if c in row and not str(row[c]).strip():
                errors.append(f"필수 열 {c}이(가) 비었습니다")
        if "sampleCollectionDate" in row and rec["meta"].get("date") and not row["sampleCollectionDate"]:
            errors.append(f"수집일 '{rec['meta'].get('date')}'을(를) YYYY / YYYY-MM / YYYY-MM-DD로 바꿀 수 없습니다")
        a = _check_authors(row.get("authors", ""))
        if a and "authors" in row:
            errors.append(a)
        if row.get("geoLocCountry") and not export.region_of(row["geoLocCountry"]):
            warns.append(f"국가 '{row['geoLocCountry']}'이(가) 알려진 이름이 아닙니다 (Pathoplexus는 국가를 표준 이름으로 받습니다)")
        if row.get("host") and not str(row["host"]).strip():
            errors.append("host가 비었습니다 (모르면 'unidentified')")
        dup = pubmd5.get(md5.get(rec["id"]))
        if dup:
            warns.append(f"같은 서열이 이미 공개 데이터에 있습니다 ({dup}) — 이미 공개된 것을 다시 제출하는 것은 아닌지 확인하세요")
        st = {t: v["status"] for t, v in rec["qc"].items()}
        if not st:
            warns.append("QC를 하지 않았습니다 (비공개 탭에서 QC 실행)")
        elif "bad" in st.values() or "error" in st.values():
            warns.append("QC에서 불량/오류가 있습니다: " + ", ".join(f"{t} {s}" for t, s in st.items() if s in ("bad", "error")))
        if rec["id"] in ex:
            warns.append("데이터셋에서 제외로 표시된 유전체입니다")
        out.append({"id": rec["id"], "row": row, "errors": errors, "warnings": warns, "status": "error" if errors else ("warn" if warns else "ok"),
                    "submitted": pkg.get(rec["id"]), "version": rec.get("version", 1), "file_path": rec["file_path"]})
    return {"organism": org, "template": tmpl, "required": [c for c in REQUIRED if c in tmpl], "rows": out}


def _terms(terms: str, until: str) -> tuple[str, str]:
    terms = (terms or "OPEN").upper()
    if terms not in ("OPEN", "RESTRICTED"):
        raise ValueError("데이터 이용 조건은 OPEN 또는 RESTRICTED입니다.")
    if terms == "OPEN":
        return "OPEN", ""
    try:
        d = date.fromisoformat((until or "").strip())
    except ValueError:
        raise ValueError("RESTRICTED는 제한이 풀리는 날짜(YYYY-MM-DD)가 필요합니다.") from None
    if d <= date.today() or d > date.today() + timedelta(days=366):
        raise ValueError("제한 기간은 오늘 이후 최대 1년까지입니다 (Pathoplexus 규칙).")
    return "RESTRICTED", d.isoformat()


def build_package(watch_id: int, ids: list[str], terms: str = "OPEN", restricted_until: str = "", group_id: str = "", note: str = "", allow_warnings: bool = True) -> dict:
    """Write the submission package of the chosen private genomes. Genomes with errors are refused (nothing is written); warnings are allowed unless asked otherwise."""
    cand = candidates(watch_id)
    if not cand["organism"]:
        raise ValueError("이 바이러스에 연결된 Pathoplexus organism을 알 수 없습니다 (지금 확인을 한 번 하면 연결됩니다).")
    terms, until = _terms(terms, restricted_until)
    chosen = [r for r in cand["rows"] if r["id"] in set(ids)]
    if not chosen:
        raise ValueError("제출할 유전체를 고르세요.")
    missing = sorted(set(ids) - {r["id"] for r in chosen})
    if missing:
        raise ValueError(f"비공개 저장소에 없는 id: {', '.join(missing[:5])}")
    bad = [r for r in chosen if r["errors"]]
    if bad:
        raise ValueError("오류가 있는 유전체는 제출 패키지에 넣을 수 없습니다:\n" + "\n".join(f"- {r['id']}: {'; '.join(r['errors'])}" for r in bad[:10]))
    if not allow_warnings and any(r["warnings"] for r in chosen):
        raise ValueError("경고가 있는 유전체가 있습니다. 확인한 뒤 다시 시도하세요.")
    con = connect()
    try:
        w = con.execute("SELECT name FROM watchlist WHERE id=?", (watch_id,)).fetchone()
    finally:
        con.close()
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    out = packages_dir() / export.slug(w["name"]) / stamp
    out.mkdir(parents=True, exist_ok=True)
    tmpl = cand["template"]
    cols = [c for c in tmpl if any(str(r["row"].get(c, "")).strip() for r in chosen) or c in REQUIRED or c == "id"]
    h_m, h_f = hashlib.sha256(), hashlib.sha256()
    with open(out / "metadata.tsv", "w", encoding="utf-8", newline="\n") as f:
        for chunk in export.iter_metadata(cols, [r["row"] for r in chosen]):
            f.write(chunk)
            h_m.update(chunk.encode("utf-8"))
    with open(out / "sequences.fasta", "w", encoding="utf-8", newline="\n") as f:
        for chunk in export.iter_private_fasta({r["id"]: r["file_path"] for r in chosen}):
            f.write(chunk)
            h_f.update(chunk.encode("utf-8"))
    sub = {"organism": cand["organism"], "created": datetime.now().isoformat(timespec="seconds"), "virus": w["name"], "genomes": len(chosen), "ids": [r["id"] for r in chosen],
           "private_versions": {r["id"]: r["version"] for r in chosen}, "data_use_terms": terms, "restricted_until": until, "group_id": group_id,
           "template_columns": tmpl, "columns_written": cols, "warnings": {r["id"]: r["warnings"] for r in chosen if r["warnings"]},
           "files": {"metadata.tsv": {"sha256": h_m.hexdigest()}, "sequences.fasta": {"sha256": h_f.hexdigest()}},
           "note": "Prepared by Virus Genome Metadata Curator. Nothing has been uploaded; upload is done by the person (README.txt)."}
    (out / "submission.json").write_text(json.dumps(sub, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (out / "submit.sh").write_text(_submit_sh(cand["organism"], terms, until, group_id), encoding="utf-8")
    (out / "submit.sh").chmod(0o750)
    (out / "README.txt").write_text(_readme(sub, cand["organism"]), encoding="utf-8")
    con = connect()
    try:
        pid = con.execute("INSERT INTO submission_package(watch_id,organism,created,dir,terms,restricted_until,group_id,n,ids_json,status,note,sha_json) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",
                          (watch_id, cand["organism"], sub["created"], str(out), terms, until, group_id, len(chosen), json.dumps(sub["ids"]), "prepared", note.strip(),
                           json.dumps({k: v["sha256"] for k, v in sub["files"].items()}))).lastrowid
        audit(con, "pathoplexus_package", f"{w['name']}: {len(chosen)}건 제출 파일 준비 ({cand['organism']}, {terms}{' ' + until if until else ''}) -> {out}")
        con.commit()
    finally:
        con.close()
    return {"id": pid, "dir": str(out), "genomes": len(chosen), "terms": terms, "restricted_until": until, "files": ["metadata.tsv", "sequences.fasta", "submission.json", "submit.sh", "README.txt"],
            "warnings": sub["warnings"]}


def _submit_sh(org: str, terms: str, until: str, group: str) -> str:
    q = f"groupId=${{GROUP_ID:-{group or 'YOUR_GROUP_ID'}}}&dataUseTermsType={terms}" + (f"&restrictedUntil={until}" if until else "")
    return f"""#!/usr/bin/env bash
# Upload this package to Pathoplexus with YOUR OWN token. The Curator never runs this file and never sees the token.
# Open data is released to INSDC and cannot be taken back; restricted data can only be shortened or opened later. Read README.txt first.
# 1) Try the demo instance first:   BACKEND=https://backend-demo.pathoplexus.org  (needs an account there)
# 2) Get a token as described at https://pathoplexus.org/docs/how-to/authentication-api and export it:  export PATHOPLEXUS_TOKEN=...
# 3) Your group id is on your group's page (reachable from your user page):  export GROUP_ID=...
set -euo pipefail
cd "$(dirname "$0")"
BACKEND="${{BACKEND:-{BACKEND}}}"
: "${{PATHOPLEXUS_TOKEN:?set PATHOPLEXUS_TOKEN (see the comments above)}}"
read -r -p "Upload to $BACKEND ({org}, {terms}{' until ' + until if until else ''})? Type 'upload' to continue: " ok
[ "$ok" = "upload" ] || {{ echo "stopped"; exit 1; }}
curl --fail-with-body -sS -X POST "$BACKEND/{org}/submit?{q}" \\
  -H "Authorization: Bearer $PATHOPLEXUS_TOKEN" \\
  -F "metadataFile=@metadata.tsv" -F "sequenceFile=@sequences.fasta"
echo
echo "The submission is now waiting for REVIEW and APPROVAL on pathoplexus.org (nothing is released before you approve)."
"""


def _readme(sub: dict, org: str) -> str:
    return f"""Pathoplexus submission package - {sub['virus']} ({org}) - {sub['genomes']} genomes
Prepared {sub['created']} by Virus Genome Metadata Curator. NOTHING HAS BEEN UPLOADED.

Data use terms: {sub['data_use_terms']}{' until ' + sub['restricted_until'] if sub['restricted_until'] else ''}
  OPEN: released without restriction and passed on to INSDC (GenBank). It cannot be taken back.
  RESTRICTED: up to one year; you can later shorten the restriction or open the data, never extend it.

Before you upload
  - You are allowed to publish these sequences (data owner / submitting group agree, patient privacy and consent are handled).
  - The authors and affiliations in metadata.tsv are the ones that should be credited.
  - Check metadata.tsv, submission.json (warnings) and sequences.fasta.

Upload on the website (recommended: it has a review step)
  1. https://pathoplexus.org - log in, choose your group, "Submit sequences", organism {org}.
  2. Select metadata.tsv and sequences.fasta, choose the same data use terms, "Upload and proceed to Approval".
  3. Review the processed result (errors are shown per sequence), fix and re-upload if needed, then approve.
  Docs: https://pathoplexus.org/docs/how-to/upload-sequences , https://pathoplexus.org/docs/how-to/approve-submissions

Upload with the API (optional)
  submit.sh does the same POST to <backend>/{org}/submit with YOUR token from the environment. Try backend-demo.pathoplexus.org first.
  A submission is still reviewed and approved afterwards on the website.

Afterwards
  In the Curator's Pathoplexus tab mark the package as submitted and, once released, record the accessions Pathoplexus gave.
  When the data are released they reach NCBI and come back through the normal 'check'; the private copy can then be left out of datasets.
Files: metadata.tsv ({sub['files']['metadata.tsv']['sha256'][:16]}...), sequences.fasta ({sub['files']['sequences.fasta']['sha256'][:16]}...)
"""


# ------------------------------------------------------------------------------------------------ the record of what was prepared / submitted
def list_packages(watch_id: int | None = None) -> list[dict]:
    con = connect()
    try:
        q = "SELECT * FROM submission_package" + (" WHERE watch_id=?" if watch_id else "") + " ORDER BY id DESC"
        out = []
        for r in con.execute(q, (watch_id,) if watch_id else ()):
            d = dict(r)
            d["ids"] = json.loads(d.pop("ids_json") or "[]")
            d["accessions"] = json.loads(d.pop("accessions_json") or "{}")
            d["sha"] = json.loads(d.pop("sha_json") or "{}")
            out.append(d)
        return out
    finally:
        con.close()


def get_package(pid: int) -> dict | None:
    return next((p for p in list_packages() if p["id"] == pid), None)


def mark(pid: int, status: str, accessions: dict | None = None, note: str | None = None) -> dict:
    """The person records what happened: 'submitted' (uploaded and approved), 'published' (with accessions), or back to 'prepared'."""
    if status not in ("prepared", "submitted", "published"):
        raise ValueError("상태는 prepared / submitted / published입니다.")
    p = get_package(pid)
    if not p:
        raise ValueError("패키지를 찾을 수 없습니다.")
    acc = {k: str(v).strip() for k, v in (accessions or {}).items() if k in p["ids"] and str(v).strip()}
    con = connect()
    try:
        con.execute("UPDATE submission_package SET status=?, accessions_json=?, note=COALESCE(?, note), submitted_at=? WHERE id=?",
                    (status, json.dumps({**p["accessions"], **acc}), note, datetime.now().isoformat(timespec="seconds") if status != "prepared" else None, pid))
        audit(con, "pathoplexus_mark", f"패키지 {pid}: {status}" + (f", 등록번호 {len(acc)}건" if acc else ""))
        con.commit()
    finally:
        con.close()
    return get_package(pid)


def zip_package(pid: int) -> bytes:
    import io
    import zipfile
    p = get_package(pid)
    if not p or not Path(p["dir"]).exists():
        raise ValueError("패키지 폴더가 없습니다.")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(Path(p["dir"]).iterdir()):
            if f.is_file():
                z.write(f, f.name)
    return buf.getvalue()
