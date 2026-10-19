"""Hand the curated genomes to other programs: a metadata table and a FASTA file that line up row for row.

Profiles (what the columns look like):
* `nextstrain` - the columns the Nextstrain mpox workflow (github.com/nextstrain/mpox) works with: `PPX_accession` as the id (its
  include/exclude lists are mapped to PPX accessions and `augur filter --exclude` compares them with the id column, so an INSDC id would
  silently skip the exclusions), standardised dates (`2023-XX-XX`), region / country / division / location, host, clade / lineage /
  outbreak and the Nextclade QC columns (`QC_rare_mutations`, ...) that its default filter reads. Genomes without a Pathoplexus accession
  cannot be identified that way and are left out (counted in the manifest). The generated config only overrides `inputs`: pass it after
  the workflow's own defaults/<build>/config.yaml (the Snakefile loads those defaults only when no config is given).
* `plain` - every curated column of the `genome` table plus the QC columns, nothing renamed.

Only open data leaves here: the curated set is defined by NCBI (INSDC) accessions, so Pathoplexus restricted-use records are never in it.
The same functions serve the CLI (`genomecurator export`), the download endpoints (`/api/export/...`) and the Export tab."""
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path

from . import __version__, completeness, config, db, log, qc, reconcile
from .db import connect

_log = log.get("export")
OPEN_TERMS_URL = "https://pathoplexus.org/about/terms-of-use/open-data"

NEXTSTRAIN_COLUMNS = [
    "PPX_accession", "PPX_accession__url", "PPX_accession_version", "PPX_accession_version__url",
    "INSDC_accession", "INSDC_accession_version", "INSDC_accession_version__url",
    "strain", "date", "region", "country", "division", "location", "length", "host", "date_released", "date_updated", "sra_accession",
    "authors", "full_authors", "institution", "dataUseTerms", "dataUseTerms__url", "restrictedUntil", "displayName", "Pathoplexus_group", "Pathoplexus_group_id", "Pathoplexus_group__url", "submission_database",
    "clade", "outbreak", "lineage", "coverage", "missing_data", "divergence", "nonACGTN", "private_substitutions",
    "QC_missing_data", "QC_mixed_sites", "QC_rare_mutations", "QC_frame_shifts", "QC_stop_codons", "frame_shifts", "is_reverse_complement"]
NS_QC = ["coverage", "missing_data", "divergence", "nonACGTN", "private_substitutions", "QC_missing_data", "QC_mixed_sites", "QC_rare_mutations",
         "QC_frame_shifts", "QC_stop_codons", "frame_shifts", "is_reverse_complement"]

_REGIONS = {
    "Africa": "algeria angola benin botswana burkina faso burundi cameroon cape verde central african republic chad comoros congo "
              "democratic republic of the congo cote d'ivoire djibouti egypt equatorial guinea eritrea eswatini ethiopia gabon gambia ghana "
              "guinea guinea-bissau kenya lesotho liberia libya madagascar malawi mali mauritania mauritius morocco mozambique namibia niger "
              "nigeria rwanda sao tome and principe senegal seychelles sierra leone somalia south africa south sudan sudan tanzania togo "
              "tunisia uganda zambia zimbabwe ivory coast zaire republic of the congo",
    "Asia": "afghanistan armenia azerbaijan bahrain bangladesh bhutan brunei cambodia china cyprus georgia india indonesia iran iraq israel japan "
            "jordan kazakhstan kuwait kyrgyzstan laos lebanon malaysia maldives mongolia myanmar nepal north korea oman pakistan palestine "
            "philippines qatar saudi arabia singapore south korea sri lanka syria taiwan tajikistan thailand timor-leste turkey turkmenistan "
            "united arab emirates uzbekistan vietnam viet nam yemen hong kong macau korea republic of korea russia",
    "Europe": "albania andorra austria belarus belgium bosnia and herzegovina bulgaria croatia czechia czech republic denmark estonia finland france "
              "germany greece hungary iceland ireland italy kosovo latvia liechtenstein lithuania luxembourg malta moldova monaco montenegro "
              "netherlands north macedonia norway poland portugal romania san marino serbia slovakia slovenia spain sweden switzerland ukraine "
              "united kingdom uk england scotland wales northern ireland vatican city russian federation",
    "North America": "antigua and barbuda bahamas barbados belize canada costa rica cuba dominica dominican republic el salvador grenada guatemala "
                     "haiti honduras jamaica mexico nicaragua panama saint kitts and nevis saint lucia saint vincent and the grenadines "
                     "trinidad and tobago usa united states puerto rico",
    "South America": "argentina bolivia brazil chile colombia ecuador guyana paraguay peru suriname uruguay venezuela",
    "Oceania": "australia fiji kiribati marshall islands micronesia nauru new zealand palau papua new guinea samoa solomon islands tonga tuvalu vanuatu",
}
# multi-word names are matched as whole strings, so the word lists above are cut at the known names instead
_NAMES = {
    "Africa": ["burkina faso", "cape verde", "central african republic", "democratic republic of the congo", "republic of the congo", "cote d'ivoire", "equatorial guinea",
               "sao tome and principe", "sierra leone", "south africa", "south sudan", "ivory coast", "guinea-bissau"],
    "Asia": ["north korea", "saudi arabia", "south korea", "sri lanka", "timor-leste", "united arab emirates", "viet nam", "hong kong", "republic of korea", "russian federation"],
    "Europe": ["bosnia and herzegovina", "czech republic", "north macedonia", "san marino", "united kingdom", "vatican city", "northern ireland"],
    "North America": ["antigua and barbuda", "costa rica", "dominican republic", "el salvador", "saint kitts and nevis", "saint lucia",
                      "saint vincent and the grenadines", "trinidad and tobago", "united states", "puerto rico"],
    "South America": [],
    "Oceania": ["marshall islands", "new zealand", "papua new guinea", "solomon islands"],
}


def _build_region_map() -> dict[str, str]:
    m = {}
    for region, text in _REGIONS.items():
        rest = text
        for n in sorted(_NAMES[region], key=len, reverse=True):
            if n in rest:
                m[n] = region
                rest = rest.replace(n, " ")
        for w in rest.split():
            m.setdefault(w, region)
    m["georgia"] = "Asia"
    m["guinea"] = "Africa"
    m["congo"] = "Africa"
    m["sudan"] = "Africa"
    m["korea"] = "Asia"
    m["russia"] = "Europe"
    return m


REGION_OF = _build_region_map()


def region_of(country: str | None) -> str:
    c = (country or "").strip().casefold()
    return REGION_OF.get(c) or REGION_OF.get(reconcile.COUNTRY_ALIAS.get(c, ""), "")


def std_date(v) -> str:
    """Nextstrain's date convention: YYYY-MM-DD, with XX for the parts that are not known ('2023-XX-XX', '2023-06-XX'); '' when unusable."""
    parts = reconcile._date_parts(str(v or "").strip())
    if not parts or not parts[0]:
        return ""
    y, m, d, rng = parts
    if rng:
        return ""
    return f"{y}-{m or 'XX'}-{d or 'XX'}"


def _iso(v) -> str:
    return str(v or "")[:10]


def _geo(r: dict) -> tuple[str, str]:
    """(division, location) from geographic_location 'USA, NY, New York City' / 'Germany' / 'Uganda: Kampala'."""
    g = str(r.get("geographic_location") or "").strip()
    if not g:
        return "", ""
    parts = [x.strip() for x in re.split(r"[,:]", g) if x.strip()]
    country = str(r.get("isolation_country") or "").strip().casefold()
    if parts and parts[0].casefold() == country:
        parts = parts[1:]
    return (parts[0] if parts else ""), (parts[1] if len(parts) > 1 else "")


# ------------------------------------------------------------------------------------------------ selection
def _load(watch_id: int) -> tuple[dict, list[dict]]:
    con = connect()
    try:
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        if not w:
            raise ValueError("바이러스를 찾을 수 없습니다.")
        rows = [dict(r) for r in con.execute(
            "SELECT g.*, s.file_path, s.md5 seq_md5, s.length seq_length FROM genome g LEFT JOIN genome_sequence s ON s.genome_id=g.genome_id "
            "WHERE g.watch_id=? AND IFNULL(g.genome_status,'')!='Removed' ORDER BY g.accession", (watch_id,))]
        qcs: dict[str, dict] = {}
        for q in con.execute("SELECT accession, tool, status, clade, lineage, metrics_json FROM qc_result WHERE watch_id=?", (watch_id,)):
            d = qcs.setdefault(q["accession"], {})
            d[q["tool"]] = {"status": q["status"], "clade": q["clade"], "lineage": q["lineage"], "m": json.loads(q["metrics_json"] or "{}")}
        ref = completeness.reference_map(con).get(watch_id)
    finally:
        con.close()
    for r in rows:
        r["_qc"] = qcs.get(r["accession"], {})
        r["completeness"], _ = completeness.classify(r.get("genome_status"), r.get("genome_length"), ref and ref[1])
    return dict(w), rows


def select(rows: list[dict], require_sequence: bool = True, min_length: int | None = None, completeness_in: list[str] | None = None,
           qc_not: list[str] | None = None, clade: str | None = None) -> list[dict]:
    out = []
    for r in rows:
        if require_sequence and not (r.get("file_path") and config.data_path(r["file_path"]).exists()):
            continue
        if min_length and (r.get("genome_length") or 0) < min_length:
            continue
        if completeness_in and r.get("completeness") not in completeness_in:
            continue
        if qc_not:
            worst = max((qc.RANK.get(v["status"], 0) for v in r["_qc"].values()), default=-1)
            names = {k for k, v in qc.RANK.items() if v == worst}
            if names & set(qc_not):
                continue
        if clade and clade not in (r.get("clade"), (r["_qc"].get("nextclade") or {}).get("clade")):
            continue
        out.append(r)
    return out


# ------------------------------------------------------------------------------------------------ rows
def _id(r: dict) -> str:
    return f"{r['accession']}.{r['accession_version']}" if r.get("accession_version") else r["accession"]


def nextstrain_row(r: dict) -> dict:
    nc = r["_qc"].get("nextclade") or {}
    ns = (nc.get("m") or {}).get("ns") or {}
    date = std_date(r.get("collection_date"))
    country = (r.get("isolation_country") or "").strip()
    division, location = _geo(r)
    ppx = str(r.get("pathoplexus_accession") or "")
    ppx_base = ppx.split(".")[0] if ppx else ""
    gid = _id(r)
    row = {
        "PPX_accession": ppx_base, "PPX_accession__url": f"https://pathoplexus.org/seq/{ppx_base}" if ppx_base else "",
        "PPX_accession_version": ppx, "PPX_accession_version__url": f"https://pathoplexus.org/seq/{ppx}" if ppx else "",
        "INSDC_accession": r["accession"], "INSDC_accession_version": gid, "INSDC_accession_version__url": f"https://www.ncbi.nlm.nih.gov/nuccore/{gid}",
        "strain": r.get("strain") or gid, "date": date, "region": region_of(country), "country": country, "division": division, "location": location,
        "length": r.get("genome_length") or "", "host": r.get("host_name") or "", "date_released": _iso(r.get("release_date")), "date_updated": _iso(r.get("update_date")),
        "sra_accession": r.get("sra_accession") or "", "authors": "?", "full_authors": "?", "institution": r.get("submitter") or "",
        "dataUseTerms": "OPEN", "dataUseTerms__url": OPEN_TERMS_URL, "restrictedUntil": "",
        "displayName": f"{(country or 'Unknown').replace(' ', '_')}/{gid}/{date or 'XXXX-XX-XX'}", "submission_database": "INSDC",
        # what the Nextstrain mpox ingest writes for records that came from INSDC (all of ours: the set is defined by NCBI)
        "Pathoplexus_group": "INSDC", "Pathoplexus_group_id": "1", "Pathoplexus_group__url": "https://pathoplexus.org/group/1",
        # curated values win; Nextclade's call is the fallback
        "clade": r.get("clade") or nc.get("clade") or "", "lineage": r.get("lineage") or nc.get("lineage") or "", "outbreak": r.get("outbreak") or ns.get("outbreak") or "",
    }
    for k in NS_QC:
        v = ns.get(k)
        row[k] = "" if v in (None, "") else v
    return row


def plain_row(r: dict) -> dict:
    row = {k: v for k, v in r.items() if not k.startswith("_") and k not in ("file_path", "seq_md5", "seq_length", "taxon_lineage_ids")}
    row["INSDC_accession_version"] = _id(r)
    st = r["_qc"].get("seqstats") or {}
    nc = r["_qc"].get("nextclade") or {}
    row.update(qc_seqstats=st.get("status", ""), qc_nextclade=nc.get("status", ""), nextclade_clade=nc.get("clade") or "", nextclade_lineage=nc.get("lineage") or "")
    return row


PROFILES = {
    "nextstrain": {"label": "Nextstrain (mpox 워크플로 형식)", "columns": NEXTSTRAIN_COLUMNS, "row": nextstrain_row, "id": "PPX_accession"},
    "plain": {"label": "일반 (genome 표 전체 + QC)", "columns": None, "row": plain_row, "id": "INSDC_accession_version"},
}


def _keep_identified(profile: str, sel: list[dict], out: list[dict]) -> tuple[list[dict], dict[str, dict], int, int]:
    """Rows need a unique, non-empty id (the profile's id column) to line up with the FASTA; the others are left out and counted."""
    idc = PROFILES[profile]["id"]
    kept, by_id, no_id, dup = [], {}, 0, 0
    for g, row in zip(sel, out):
        k = str(row.get(idc) or "")
        if not k:
            no_id += 1
        elif k in by_id:
            dup += 1
        else:
            kept.append(row)
            by_id[k] = g
    return kept, by_id, no_id, dup


def table(watch_id: int, profile: str = "nextstrain", **filters) -> tuple[list[str], list[dict], dict, dict]:
    """(columns, rows, selected genome rows keyed by the id column, watch)."""
    if profile not in PROFILES:
        raise ValueError(f"알 수 없는 형식입니다: {profile} (가능: {', '.join(PROFILES)})")
    w, rows = _load(watch_id)
    sel = select(rows, **filters)
    prof = PROFILES[profile]
    out = [prof["row"](r) for r in sel]
    out, by_id, _, _ = _keep_identified(profile, sel, out)
    cols = prof["columns"] or (list(out[0]) if out else ["INSDC_accession_version"])
    return cols, out, by_id, w


def preview(watch_id: int, profile: str = "nextstrain", **filters) -> dict:
    """Counts and warnings, without writing anything."""
    if profile not in PROFILES:
        raise ValueError(f"알 수 없는 형식입니다: {profile} (가능: {', '.join(PROFILES)})")
    w, rows = _load(watch_id)
    sel = select(rows, **filters)
    out_all = [PROFILES[profile]["row"](r) for r in sel]
    out, by_id, no_id, dup = _keep_identified(profile, sel, out_all)
    n = len(out)
    warn = []
    no_seq = sum(1 for r in rows if not (r.get("file_path") and config.data_path(r["file_path"]).exists()))
    if filters.get("require_sequence", True) and no_seq:
        warn.append(f"서열 파일이 없는 유전체 {no_seq:,}건은 제외됩니다 (유전체 화면의 Get FASTA로 받으면 포함됩니다).")
    idc = PROFILES[profile]["id"]
    if no_id:
        warn.append(f"`{idc}` 값이 없는 {no_id:,}건은 제외됩니다" + (" (Nextstrain mpox 워크플로는 Pathoplexus accession을 id로 쓰고 포함·제외 목록도 그것으로 맞춥니다)." if profile == "nextstrain" else "."))
    if dup:
        warn.append(f"`{idc}`가 겹치는 {dup:,}건은 제외됩니다 (id는 유일해야 합니다).")
    if profile == "nextstrain" and n:
        miss_qc = sum(1 for r in out if not r["QC_rare_mutations"])
        if miss_qc:
            warn.append(f"Nextclade QC 값이 없는 {miss_qc:,}건: Nextstrain mpox의 기본 필터가 `QC_rare_mutations`를 보므로 이 행들은 워크플로에서 걸러집니다. QC 탭에서 QC를 (다시) 실행하세요.")
        for col, label in (("date", "수집일"), ("region", "대륙(region)"), ("country", "국가")):
            m = sum(1 for r in out if not r[col])
            if m:
                warn.append(f"{label}이(가) 비었거나 형식이 맞지 않는 {m:,}건 (Nextstrain에서는 빈 값으로 처리됩니다).")
    return {"virus": w["name"], "profile": profile, "genomes": len(rows), "selected": n, "left_out_no_id": no_id, "left_out_duplicate": dup, "warnings": warn}


# ------------------------------------------------------------------------------------------------ writing
def _tsv_cell(v) -> str:
    return str(v if v is not None else "").replace("\t", " ").replace("\n", " ").replace("\r", " ")


def iter_metadata(cols: list[str], rows: list[dict]):
    yield "\t".join(cols) + "\n"
    for r in rows:
        yield "\t".join(_tsv_cell(r.get(c)) for c in cols) + "\n"


def _seq(path: Path) -> str:
    seq, started = [], False
    with open(path, encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith(">"):
                if started:
                    break
                started = True
            else:
                seq.append(ln.strip())
    return "".join(seq)


def iter_fasta(profile: str, sel: dict[str, dict], width: int = 70):
    """One record per selected genome, id = the profile's id column; sequence read from the stored file (header replaced)."""
    for gid, r in sel.items():
        if not (r.get("file_path") and config.data_path(r["file_path"]).exists()):
            continue
        s = _seq(config.data_path(r["file_path"]))
        yield f">{gid}\n" + "\n".join(s[i:i + width] for i in range(0, len(s), width)) + "\n"


def slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "virus"


def exports_dir() -> Path:
    return config.DATA_DIR / "exports"


def write(watch_id: int, out_dir: Path | None = None, profile: str = "nextstrain", progress=None, nextstrain_dir: Path | None = None, **filters) -> dict:
    """Write metadata.tsv, sequences.fasta, manifest.json and nextstrain-config.yaml. Default folder: data/exports/<virus>/<timestamp>/.
    With nextstrain_dir (the `phylogenetic/` folder of a nextstrain/mpox checkout) the files go to <it>/data/genomecurator/ and the config to
    <it>/build-configs/genomecurator/config.yaml with paths relative to that folder - the form that also works inside the Docker runtime."""
    p = progress or (lambda *a: None)
    p(5, "선택한 유전체와 QC 결과를 모으는 중…")
    cols, rows, sel, w = table(watch_id, profile, **filters)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    ns = Path(nextstrain_dir).expanduser() if nextstrain_dir else None
    if ns and not (ns / "Snakefile").exists():
        raise ValueError(f"{ns}에 Snakefile이 없습니다. nextstrain/mpox 저장소의 phylogenetic/ 폴더를 지정하세요.")
    out = (ns / "data" / "genomecurator") if ns else (Path(out_dir) if out_dir else exports_dir() / slug(w["name"]) / stamp)
    out.mkdir(parents=True, exist_ok=True)
    meta, fa = out / "metadata.tsv", out / "sequences.fasta"
    h_meta, h_fa = hashlib.sha256(), hashlib.sha256()
    with open(meta, "w", encoding="utf-8", newline="\n") as f:
        for chunk in iter_metadata(cols, rows):
            f.write(chunk)
            h_meta.update(chunk.encode("utf-8"))
    p(40, f"서열 {len(sel):,}건 쓰는 중…")
    n_seq = 0
    with open(fa, "w", encoding="utf-8", newline="\n") as f:
        for chunk in iter_fasta(profile, sel):
            f.write(chunk)
            h_fa.update(chunk.encode("utf-8"))
            n_seq += 1
            if n_seq % 500 == 0:
                p(40 + 55 * n_seq / max(1, len(sel)), f"서열 {n_seq:,}/{len(sel):,}")
    pv = preview(watch_id, profile, **filters)
    con = connect()
    try:
        tags = sorted({f"{r[0]}@{r[1]}" for r in con.execute("SELECT DISTINCT dataset, dataset_tag FROM qc_result WHERE watch_id=? AND tool='nextclade'", (watch_id,))})
        nver = sorted({r[0] for r in con.execute("SELECT DISTINCT version FROM qc_result WHERE watch_id=? AND tool='nextclade'", (watch_id,)) if r[0]})
        dbv = db.version_info()
    finally:
        con.close()
    manifest = {"generated": datetime.now().isoformat(timespec="seconds"), "virus": w["name"], "taxon": w["taxon"], "profile": profile,
                "app": {"name": "genomecurator", "version": __version__, "db_schema": dbv.get("db_schema_version")},
                "filters": {k: v for k, v in filters.items() if v not in (None, [], "")}, "rows": len(rows), "sequences": n_seq,
                "id_column": PROFILES[profile]["id"], "warnings": pv["warnings"], "nextclade": {"datasets": tags, "versions": nver},
                "files": {"metadata.tsv": {"sha256": h_meta.hexdigest(), "columns": cols}, "sequences.fasta": {"sha256": h_fa.hexdigest()}},
                "data_use": "Open data only: the curated set is defined by INSDC (NCBI) accessions, so Pathoplexus restricted-use records are not included. "
                            "Cite the original submitters and sources (NCBI, Pathoplexus, BV-BRC); see docs/LICENSES.md."}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if ns:
        cfg_path = ns / "build-configs" / "genomecurator" / "config.yaml"
        cfg_path.parent.mkdir(parents=True, exist_ok=True)
        cfg = nextstrain_config(out, profile, relative_to=ns)
    else:
        cfg_path = out / "nextstrain-config.yaml"
        cfg = nextstrain_config(out, profile)
    cfg_path.write_text(cfg, encoding="utf-8")
    con = connect()
    try:
        db.audit(con, "export", f"{w['name']}: {len(rows)}행, 서열 {n_seq}건 -> {out} ({profile})")
        con.commit()
    finally:
        con.close()
    p(100, "내보내기 완료")
    _log.info("export %r profile=%s rows=%d sequences=%d -> %s", w["name"], profile, len(rows), n_seq, out)
    return {"dir": str(out), "rows": len(rows), "sequences": n_seq, "warnings": pv["warnings"], "config": str(cfg_path),
            "files": ["metadata.tsv", "sequences.fasta", "manifest.json"] + ([] if ns else ["nextstrain-config.yaml"])}


def nextstrain_config(out: Path, profile: str = "nextstrain", relative_to: Path | None = None) -> str:
    """A config *overlay* for the Nextstrain mpox workflow: only `inputs` (local files instead of data.nextstrain.org) and use_remote_files.
    The workflow's Snakefile loads its own defaults/<build>/config.yaml only when no config is given, so pass the build's defaults first."""
    def path(name: str) -> str:
        f = out / name
        return (f.relative_to(relative_to) if relative_to else f).as_posix()
    cfgp = "build-configs/genomecurator/config.yaml" if relative_to else str(out / "nextstrain-config.yaml")
    idc = PROFILES[profile]["id"]
    return (f"# Overlay for github.com/nextstrain/mpox (phylogenetic/). Run with the build's own defaults FIRST, then this file:\n"
            f"#   cd phylogenetic\n"
            f"#   nextstrain build . --configfile defaults/mpxv/config.yaml {cfgp}        # all clades\n"
            f"#   nextstrain build . --configfile defaults/hmpxv1/config.yaml {cfgp}      # clade IIb   (or defaults/clade-i/config.yaml)\n"
            f"# (the Snakefile only loads defaults/*/config.yaml when no config is given, so this file alone is not a complete config)\n"
            f"# Inputs replace the default URLs; local paths are used as they are, nothing is downloaded.\n"
            f"{'' if relative_to else '# Absolute paths do not exist inside the Docker runtime: use --nextstrain-dir to write into the checkout, or run with --ambient/--conda.' + chr(10)}"
            f"use_remote_files: false\n"
            f"inputs:\n  - name: genomecurator\n    metadata: \"{path('metadata.tsv')}\"\n    sequences: \"{path('sequences.fasta')}\"\n"
            + (f"strain_id_field: \"{idc}\"\n" if idc != "PPX_accession" else ""))


def list_exports(watch_name: str | None = None, limit: int = 20) -> list[dict]:
    base = exports_dir()
    out = []
    if not base.exists():
        return out
    for vdir in sorted(base.iterdir()):
        if watch_name and vdir.name != slug(watch_name):
            continue
        for d in sorted((x for x in vdir.iterdir() if x.is_dir()), reverse=True)[:limit]:
            mf = d / "manifest.json"
            try:
                m = json.loads(mf.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            out.append({"virus": m.get("virus"), "dir": str(d), "generated": m.get("generated"), "rows": m.get("rows"), "sequences": m.get("sequences"), "profile": m.get("profile")})
    return out
