from genomecurator import curate, sources


def test_split_version():
    assert sources.split_version("PZ173267.2") == ("PZ173267", 2)
    assert sources.split_version("PZ173267") == ("PZ173267", None)


def test_merge_priority_and_provenance():
    per = {
        "ncbi": {"A1": {"host_name": "bat", "strain": "s1"}},
        "bvbrc": {"A1": {"host_name": "Bat?", "genotype": "B"}},
        "pathoplexus": {"A1": {"host_name": "Pteropus"}},
    }
    rec, prov = curate.merge(per)["A1"]
    assert rec["host_name"] == "Pteropus" and prov["host_name"] == "pathoplexus"
    assert rec["strain"] == "s1" and prov["strain"] == "ncbi"
    assert rec["genotype"] == "B" and prov["genotype"] == "bvbrc"


def test_merge_only_ncbi_accessions():
    per = {"ncbi": {"A1": {}}, "bvbrc": {"B9": {"strain": "x"}}}
    assert set(curate.merge(per)) == {"A1"}


def test_diff_detects_version_and_metadata():
    old = {"accession_version": 1, "host_name": "bat"}
    new = {"accession_version": 2, "host_name": "bat", "strain": "s"}
    d = curate._diff(old, new)
    assert d["accession_version"] == [1, 2] and "strain" in d and "host_name" not in d
    assert curate._diff(old, dict(old)) == {}


def test_merge_coerces_list_values():
    per = {"ncbi": {"A1": {}}, "bvbrc": {"A1": {"segment": ["L"], "strain": "x"}}}
    rec, _ = curate.merge(per)["A1"]
    assert rec["segment"] == "L" and rec["strain"] == "x"


def test_base_template_script_has_no_broken_string_literals():
    """A raw newline inside a "..." JS string once broke the whole script (runJob undefined => buttons dead)."""
    import re
    from pathlib import Path
    t = (Path(curate.__file__).parent / "templates" / "base.html").read_text(encoding="utf-8")
    js = re.findall(r"<script>(.*?)</script>", t, re.S)[-1]
    for n, line in enumerate(js.splitlines(), 1):
        assert len(re.findall(r'(?<!\\)"', line)) % 2 == 0, f"unbalanced quotes in base.html script line {n}: {line[:80]}"
    for name in ("runJob", "pollJob", "deploy"):
        assert f"function {name}(" in js


def test_cross_platform_helpers(tmp_path):
    from genomecurator import config
    assert config.safe_name("Monkeypox virus") == "Monkeypox_virus"
    assert config.safe_name("a/b:c*d") == "a_b_c_d"          # characters illegal on Windows
    assert config.safe_name("  ..  ") == "virus"             # never empty / never a dot-only name
    assert config.DATA_DIR.is_absolute()                     # independent of the current directory
    f = tmp_path / "x.fna"
    curate._write_fasta(f, ">AB123.1 test virus", "ACGT" * 40)
    raw = f.read_bytes()
    assert b"\r" not in raw                                  # LF on every OS (plain write_text gives CRLF on Windows)
    assert raw.splitlines()[0] == b">AB123.1 test virus" and len(raw.splitlines()[1]) == 70
    assert curate._read_staged(f, "AB123.1") == "ACGT" * 40


def test_sources_have_no_bom_and_are_utf8():
    import codecs
    from pathlib import Path
    root = Path(curate.__file__).parent
    for p in list(root.glob("*.py")) + list(root.glob("templates/*.html")):
        raw = p.read_bytes()
        assert not raw.startswith(codecs.BOM_UTF8), f"BOM in {p.name}"
        raw.decode("utf-8")


def test_stored_fasta_paths_survive_moving_the_data_folder(tmp_path, monkeypatch):
    from genomecurator import config
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    f = tmp_path / "genomes" / "Some_virus" / "1.1.fna"
    f.parent.mkdir(parents=True)
    f.write_text(">x\nACGT\n", encoding="utf-8")
    assert config.rel_data_path(f) == "genomes/Some_virus/1.1.fna"               # portable form
    assert config.data_path("genomes/Some_virus/1.1.fna") == f
    # absolute paths written by older versions (here: a Windows path on any OS) are re-rooted
    assert config.data_path(r"C:\Users\me\proj\data\genomes\Some_virus\1.1.fna") == f
    assert config.data_path("/home/me/proj/data/genomes/Some_virus/1.1.fna") == f


def test_watch_update_and_delete(tmp_path, monkeypatch):
    from genomecurator import config
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    con = connect()
    for i, n in ((1, "Alpha virus"), (2, "Beta virus")):
        con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(?,?,?,?)", (i, n, n, 100 + i))
    fa = tmp_path / "genomes" / "Alpha_virus" / "101.1.fna"
    fa.parent.mkdir(parents=True)
    fa.write_text(">x\nACGT\n", encoding="utf-8")
    con.execute("INSERT INTO genome(genome_id,watch_id,accession) VALUES('101.1',1,'AB1')")
    con.execute("INSERT INTO genome_sequence(accession,genome_id,file_path) VALUES('AB1','101.1',?)", (config.rel_data_path(fa),))
    con.execute("INSERT INTO field_provenance(genome_id,field,source) VALUES('101.1','strain','ncbi')")
    con.execute("INSERT INTO changeset(id,watch_id) VALUES(7,1)")
    con.execute("INSERT INTO changeset_item(changeset_id,accession,kind) VALUES(7,'AB1','new')")
    con.execute("INSERT INTO genome(genome_id,watch_id,accession) VALUES('102.1',2,'CD2')")
    con.commit()
    con.close()

    assert curate.update_watch(1, name="Alpha renamed", enabled=False)
    row = connect().execute("SELECT name, enabled FROM watchlist WHERE id=1").fetchone()
    assert (row["name"], row["enabled"]) == ("Alpha renamed", 0)
    import pytest
    with pytest.raises(ValueError):
        curate.update_watch(1, name="Beta virus")                    # name already used by another watch
    assert curate.watch_impact(1) == {"name": "Alpha renamed", "genomes": 1, "fasta": 1, "changesets": 1, "items": 1}

    res = curate.delete_watch(1)
    assert res["genomes"] == 1 and res["files_removed"] == 1
    assert not fa.exists() and not fa.parent.exists()                   # file and its empty folder are gone
    db = connect()
    for t in ("genome", "genome_sequence", "field_provenance", "changeset", "changeset_item"):
        n = db.execute(f"SELECT COUNT(*) FROM {t} WHERE " + {"genome": "watch_id=1", "genome_sequence": "genome_id='101.1'", "field_provenance": "genome_id='101.1'",
                                                          "changeset": "watch_id=1", "changeset_item": "changeset_id=7"}[t]).fetchone()[0]
        assert n == 0, t
    assert db.execute("SELECT COUNT(*) FROM genome WHERE watch_id=2").fetchone()[0] == 1   # other watches untouched
    assert db.execute("SELECT COUNT(*) FROM watchlist").fetchone()[0] == 1


def test_pathoplexus_slugs_are_discovered_from_the_site(monkeypatch):
    from genomecurator import sources

    class R:
        status_code = 200
        text = '<a href="/ebola-bdbv/search">x</a> <a href="/andv/search">y</a> <a href="/about/terms">z</a>'

    monkeypatch.setattr(sources, "_get", lambda c, url, **kw: R())
    slugs = sources.pathoplexus_organisms()
    assert "ebola-bdbv" in slugs and "andv" in slugs and "mpox" in slugs   # site links + built-in fallback list
    assert "about" not in slugs


def test_backup_roundtrip(tmp_path, monkeypatch):
    from genomecurator import backup, config
    from genomecurator.db import connect
    src = tmp_path / "src"
    monkeypatch.setattr(config, "DATA_DIR", src)
    monkeypatch.setattr(config, "DB_PATH", src / "g.db")
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'V','V',9)")
    con.execute("INSERT INTO genome(genome_id,watch_id,accession,assembly_method,clade) VALUES('9.1',1,'AB1','spades','Ia')")
    con.execute("INSERT INTO changeset(id,watch_id) VALUES(5,1)")
    con.execute("INSERT INTO changeset_item(changeset_id,accession,kind,decision) VALUES(5,'AB1','new','rejected'),(5,'AB2','new','pending')")
    con.execute("INSERT INTO source_record(source,accession,watch_id,payload_json) VALUES('ncbi','AB1',1,'{}')")
    fa = src / "genomes" / "V" / "9.1.fna"
    fa.parent.mkdir(parents=True)
    fa.write_text(">AB1.1\nACGT\n", encoding="utf-8")
    con.execute("INSERT INTO genome_sequence(accession,genome_id,file_path) VALUES('AB1','9.1','genomes/V/9.1.fna')")
    con.commit()
    con.close()
    info = backup.export_backup(tmp_path / "b.zip")
    assert info["fasta_files"] == 1 and info["tables"]["genome"] == 1
    dst = tmp_path / "dst"
    monkeypatch.setattr(config, "DATA_DIR", dst)
    monkeypatch.setattr(config, "DB_PATH", dst / "g.db")
    backup.restore_backup(tmp_path / "b.zip")
    db = connect()
    assert db.execute("SELECT COUNT(*) FROM genome").fetchone()[0] == 1
    assert tuple(db.execute("SELECT assembly_method, clade FROM genome").fetchone()) == ("spades", "Ia")   # migrated columns survive
    assert [r[0] for r in db.execute("SELECT decision FROM changeset_item")] == ["rejected"]   # decisions kept, pending dropped
    assert db.execute("SELECT COUNT(*) FROM source_record").fetchone()[0] == 0                   # re-downloadable data is not backed up
    assert config.data_path("genomes/V/9.1.fna").read_text(encoding="utf-8") == ">AB1.1\nACGT\n"
    import pytest
    with pytest.raises(FileExistsError):
        backup.restore_backup(tmp_path / "b.zip")


def test_delete_changeset(tmp_path, monkeypatch):
    import pytest
    from genomecurator import config, curate
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'V','V',9)")
    con.execute("INSERT INTO genome(genome_id,watch_id,accession) VALUES('9.1',1,'AB1')")
    con.execute("INSERT INTO changeset(id,watch_id,status) VALUES(5,1,'APPLIED'),(6,1,'READY')")
    con.execute("INSERT INTO changeset_item(changeset_id,accession,kind,decision) VALUES(5,'AB1','new','applied'),(5,'AB2','new','rejected'),(6,'AB3','new','pending')")
    con.commit()
    con.close()
    imp = curate.changeset_impact(5)
    assert imp["items"] == 2 and imp["rejected"] == 1
    curate.delete_changeset(5)
    db = connect()
    assert [r[0] for r in db.execute("SELECT id FROM changeset")] == [6]
    assert [r[0] for r in db.execute("SELECT changeset_id FROM changeset_item")] == [6]   # other changeset untouched
    assert db.execute("SELECT COUNT(*) FROM genome").fetchone()[0] == 1                    # applied genomes stay
    with pytest.raises(ValueError):
        curate.delete_changeset(5)


def test_reconcile_classify():
    from genomecurator.reconcile import classify
    assert classify("date", "2022-07", "2022-07-11") == "cosmetic"
    assert classify("date", "2024", "2024-01-01/2024-12-09") == "cosmetic"
    assert classify("date", "2026-03", "Mar-2026") == "cosmetic"
    assert classify("date", "2022-05-30", "1922-05-30") == "conflict"
    assert classify("place", "USA: California", "USA, California") == "cosmetic"
    assert classify("country", "Zaire", "Democratic Republic of the Congo") == "cosmetic"
    assert classify("country", "Democratic Republic of the Congo", "South Korea") == "conflict"
    assert classify("host", "Homo sapiens", "Homo sapiens; 32; male") == "cosmetic"
    assert classify("strain", "843593", "hMpxV/USA/ID-IBL-843593/2022") == "cosmetic"
    assert classify("strain", "Monkeypox_PT1009_2027", "Monkeypox_PT1009_2026") == "conflict"
    assert classify("num", 2, "1") == "conflict" and classify("num", 5, "5") == "match"
    assert classify("set", "SRR1; SRR2", "SRR2,SRR1") == "cosmetic"


def test_reconcile_analyze_and_store(tmp_path, monkeypatch):
    from genomecurator import config, reconcile
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    ps = {"ncbi": {"A1": {"genome_length": 100, "collection_date": "2022-07"}, "A2": {"genome_length": 5}},
          "bvbrc": {"A1": {"genome_length": 99, "collection_date": "2022-07-11"}, "B9": {"genome_length": 7}},
          "pathoplexus": {"A1": {"accession_version": 1}, "A3": {}}}
    rep, rows = reconcile.analyze(ps)
    assert rep["coverage"]["not_in_ncbi"] == {"bvbrc": 1, "pathoplexus": 1}
    assert rep["coverage"]["all_sources"] == 1
    assert rep["fields"]["genome_length"]["ncbi|bvbrc"]["conflict"] == 1
    assert rep["fields"]["collection_date"]["ncbi|bvbrc"]["cosmetic"] == 1
    con = connect()
    reconcile.store(con, 1, rep, rows)
    con.commit()
    con.close()
    assert reconcile.load(1)["totals"] == rep["totals"]
    assert [r["field"] for r in reconcile.conflicts(1)] == ["genome_length"]


def test_verify_detects_problems(tmp_path, monkeypatch):
    from genomecurator import config, verify
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'V','V',9)")
    con.execute("INSERT INTO genome(genome_id,watch_id,accession,accession_version,genome_length,collection_date,genome_status) VALUES('1',1,'AB1',1,8,'2999-01-01','Complete')")
    fa = tmp_path / "genomes" / "V" / "1.fna"
    fa.parent.mkdir(parents=True)
    fa.write_text(">AB1.1\nACGT\n", encoding="utf-8")
    con.execute("INSERT INTO genome_sequence(accession,genome_id,length,md5,file_path) VALUES('AB1','1',8,'x','genomes/V/1.fna')")
    con.commit()
    con.close()
    res = {c["name"]: c for c in verify.run()}
    assert res["수집일이 미래"]["count"] == 1
    assert res["FASTA 실제 길이 ≠ 기록 길이"]["count"] == 1 and res["FASTA md5 ≠ 기록 md5"]["count"] == 1
    assert res["watch가 없는 유전체"]["count"] == 0


def test_source_agreement_marks_conflicting_items(tmp_path, monkeypatch):
    from genomecurator import config, curate, reconcile
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    con = connect()
    rep, rows = reconcile.analyze({"ncbi": {"A1": {"genome_length": 1, "collection_date": "2022"}, "A2": {"genome_length": 5}},
                                   "bvbrc": {"A1": {"genome_length": 2, "collection_date": "2022-01"}, "A2": {"genome_length": 5}}})
    reconcile.store(con, 1, rep, rows)
    con.commit()
    con.close()
    ag = curate.source_agreement(1)
    assert ag["A1"]["conflict"] == {"genome_length"} and ag["A1"]["cosmetic"] == {"collection_date"}
    assert "A2" not in ag


def test_seqcheck_compare():
    from genomecurator.seqcheck import _norm, compare
    assert _norm(" acg\nT ") == "ACGT"
    assert compare("ACGTACGT", "ACGTACGT")["status"] == "same"
    d = compare("ACGTACGTAA", "ACGTTTTTAA")
    assert d["status"] == "diff" and d["prefix"] == 4 and d["suffix"] == 3
    assert compare("AAAC", "GTTT")["status"] == "diff" and "역상보" in compare("AAAC", "GTTT")["detail"]
    assert "길이" in compare("ACGT", "ACG")["detail"]


def test_seq_verdict_resolves_length_conflict(tmp_path, monkeypatch):
    from genomecurator import config, curate, reconcile
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    con = connect()
    rep, rows = reconcile.analyze({"ncbi": {"A1": {"genome_length": 1, "strain": "x1"}, "A2": {"genome_length": 1}},
                                   "bvbrc": {"A1": {"genome_length": 2, "strain": "y9"}, "A2": {"genome_length": 2}}})
    reconcile.store(con, 1, rep, rows)
    con.executemany("INSERT INTO sequence_check(watch_id,accession,source,status) VALUES(1,?,?,?)",
                    [("A1", "bvbrc", "same"), ("A2", "bvbrc", "diff")])
    con.commit()
    con.close()
    ag = curate.source_agreement(1)
    assert ag["A1"]["conflict"] == {"strain"} and ag["A1"]["resolved"] == {"genome_length"} and ag["A1"]["seq"] == "same"
    assert ag["A2"]["conflict"] == {"genome_length"} and ag["A2"]["seq"] == "diff"


def test_source_failure_stops_check_without_changes(tmp_path, monkeypatch):
    import httpx
    import pytest
    from genomecurator import config, curate, sources
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'V','V',9)")
    con.commit()
    con.close()
    monkeypatch.setattr(sources, "ncbi_species_scope", lambda tid: None)
    monkeypatch.setattr(sources, "fetch_ncbi", lambda taxon, filters=None: {"A1": {"genome_length": 5}})
    def boom(*a, **k):
        raise httpx.ConnectError("down")
    monkeypatch.setattr(sources, "fetch_pathoplexus", boom)
    monkeypatch.setattr(sources, "fetch_bvbrc", lambda tid: {})
    with pytest.raises(sources.SourceUnavailable):
        curate.check("V")
    db = connect()
    assert db.execute("SELECT COUNT(*) FROM changeset").fetchone()[0] == 0   # nothing was merged or recorded as changed
    assert db.execute("SELECT COUNT(*) FROM genome").fetchone()[0] == 0


def test_logging_context_masking_and_audit(tmp_path, monkeypatch):
    import logging
    from genomecurator import config, log, sources
    from genomecurator.db import audit_now, connect
    assert "api_key=***" in sources._mask({"api_key": "SECRET", "x": 1}) and "SECRET" not in sources._mask({"api_key": "SECRET"})
    records = []
    h = logging.Handler()
    h.emit = lambda r: records.append((r, getattr(r, "ctx", None)))
    h.addFilter(log._Context())
    lg = logging.getLogger("genomecurator")
    lg.addHandler(h)
    lg.setLevel(logging.DEBUG)
    with log.context(job="j1", watch="V"):
        log.get("t").info("hello")
    assert records[-1][1] == " [job=j1 watch=V]"
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    audit_now("x_failed", "d" * 1000)
    detail = connect().execute("SELECT detail FROM audit_log WHERE action=?", ("x_failed",)).fetchone()[0]
    assert len(detail) < 700 and "(+400자)" in detail
    lg.removeHandler(h)


def test_schema_check_detects_drift_and_keeps_baseline(tmp_path, monkeypatch):
    from genomecurator import schema_check as sc
    rows = [{"accession": "A1.1", "length": 10, "virus": {"tax_id": 9, "organism_name": "V", "lineage": [{"tax_id": 1, "name": "x"}]},
             "host": {"organism_name": "H"}, "isolate": {"collection_date": "2020"}} for _ in range(5)]
    seen, produced = sc._read_set("ncbi", rows)
    assert {"accession", "length", "host.organism_name", "virus.lineage[].tax_id"} <= seen      # the real adapter's reads are measured
    live = sc.profile(rows)
    drifted = [{k: v for k, v in r.items() if k != "host"} | {"seqLength": 10} for r in rows]
    prev = {p: dict(v) for p, v in live.items()}
    d_live = sc.profile(drifted)
    f = sc._compare("ncbi", d_live, seen, prev, {"top": ["reports"]}, ["reports"])
    kinds = {(x["kind"], x["path"]) for x in f}
    assert ("removed", "host") in kinds or ("removed", "host.organism_name") in kinds
    assert any(x["severity"] == "fail" for x in f)
    assert sc._compare("ncbi", live, seen, prev, {"top": ["reports"]}, ["reports"]) == []       # identical shape: nothing to report


def test_schema_local_check_finds_missing_column(tmp_path, monkeypatch):
    import sqlite3
    from genomecurator import config, schema_check as sc
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    connect().close()
    raw = sqlite3.connect(config.DB_PATH)
    raw.execute("ALTER TABLE genome DROP COLUMN segment")
    raw.commit()
    raw.close()
    out = sc.local_checks({})
    assert any(f["severity"] == "fail" and f["kind"] == "column" and "segment" in f["detail"] for f in out)


def test_schema_daily_modes_avoid_false_alarms():
    from genomecurator import schema_check as sc
    full = [{"accession": "A.1", "length": 5, "host": {"organism_name": "H"}} for _ in range(10)]
    seen, _ = sc._read_set("ncbi", full)
    base = sc.profile(full)
    old3 = [{"accession": "A.1", "length": 5} for _ in range(3)]          # three old records of one virus that simply have no host
    sh = {"top": ["reports"]}
    assert sc._compare("ncbi", sc.profile(old3), seen, base, sh, ["reports"], mode="virus", n=3) == []   # one virus' 3 records prove nothing
    assert sc._compare("ncbi", sc.profile(old3 + full[:3]), seen, base, sh, ["reports"], mode="pooled", n=6) == []   # present elsewhere: fine
    gone = [{"accession": "A.1", "length": 5} for _ in range(6)]
    f = sc._compare("ncbi", sc.profile(gone), seen, base, sh, ["reports"], mode="pooled", n=6)
    assert any(x["kind"] == "removed" and x["path"].startswith("host") for x in f)   # gone from every virus: real
    assert sc._compare("ncbi", sc.profile(old3), seen, base, {"top": ["results"]}, ["reports"], mode="virus", n=3)[0]["kind"] == "shape"


def test_version_is_recorded_in_db_backup_and_changeset(tmp_path, monkeypatch):
    import json
    import zipfile
    import pytest
    from genomecurator import DB_SCHEMA_VERSION, __version__, backup, config, db
    from genomecurator.db import connect, parse_version, version_info
    monkeypatch.setattr(config, "DATA_DIR", tmp_path / "a")
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "a" / "g.db")
    db._recorded.discard(str(config.DB_PATH))
    connect().close()
    meta = {r[0]: r[1] for r in connect().execute("SELECT key, value FROM app_meta")}
    assert meta["app_version"] == __version__ and meta["first_version"] == __version__ and int(meta["db_schema_version"]) == DB_SCHEMA_VERSION
    assert sum(1 for r in connect().execute("SELECT 1 FROM audit_log WHERE action='version'")) == 1      # recorded once, not on every connect
    assert version_info()["version"] == __version__
    # a database last handled by an older version: the change is written to the audit log and first_version is kept
    con = connect()
    con.execute("UPDATE app_meta SET value='0.0.9' WHERE key='app_version'")
    con.execute("UPDATE app_meta SET value='0.0.9' WHERE key='first_version'")
    con.commit()
    con.close()
    db._recorded.discard(str(config.DB_PATH))
    con = connect()
    assert [r[0] for r in con.execute("SELECT value FROM app_meta WHERE key IN ('app_version','first_version') ORDER BY key")] == [__version__, "0.0.9"]
    assert "이전 0.0.9" in con.execute("SELECT detail FROM audit_log WHERE action='version' ORDER BY id DESC").fetchone()[0]
    con.close()
    # backups carry the version; one from a newer database layout is refused
    zpath = tmp_path / "b.zip"
    info = backup.export_backup(zpath)
    assert info["app_version"] == __version__ and info["db_schema_version"] == DB_SCHEMA_VERSION
    newer = tmp_path / "newer.zip"
    with zipfile.ZipFile(zpath) as zi, zipfile.ZipFile(newer, "w") as zo:
        for n in zi.namelist():
            zo.writestr(n, json.dumps({**json.loads(zi.read(n)), "db_schema_version": DB_SCHEMA_VERSION + 1}) if n == "manifest.json" else zi.read(n))
    monkeypatch.setattr(config, "DATA_DIR", tmp_path / "b")
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "b" / "g.db")
    with pytest.raises(ValueError, match="더 새로운"):
        backup.restore_backup(newer)
    assert parse_version("0.10.2") > parse_version("0.9.9") and parse_version("v1.0") == (1, 0, 0)


def test_changelog_matches_the_code_version():
    """The release rules in docs/RELEASING.md, enforced: pyproject reads the version from the package, both changelogs (English and
    Korean) have an entry for it with a date, keep an [Unreleased] section on top, list releases newest first, and agree with each other."""
    import re
    from pathlib import Path
    from genomecurator import __version__
    root = Path(__file__).resolve().parents[1]
    per_file = {}
    for name in ("CHANGELOG.md", "CHANGELOG.ko.md"):
        log = (root / name).read_text(encoding="utf-8")
        heads = re.findall(r"^## \[([^\]]+)\](?: - (\d{4}-\d{2}-\d{2}))?", log, re.M)
        assert heads and heads[0][0] == "Unreleased", f"{name} must start with an [Unreleased] section"
        released = [h for h in heads if h[0] != "Unreleased"]
        assert released and released[0][0] == __version__, f"{name}: newest release {released[0][0] if released else None} != code version {__version__}"
        assert all(d for _, d in released), f"{name}: every released version needs a date"
        vs = [tuple(int(x) for x in v.split(".")) for v, _ in released]
        assert vs == sorted(vs, reverse=True) and len(set(vs)) == len(vs), f"{name}: versions must be unique and newest first"
        assert f"[{__version__}]:" in log, f"{name}: link reference for the current version is missing at the bottom"
        per_file[name] = released
    assert per_file["CHANGELOG.md"] == per_file["CHANGELOG.ko.md"], "English and Korean changelogs list different versions/dates"
    assert 'dynamic = ["version"]' in (root / "pyproject.toml").read_text(encoding="utf-8")


def test_every_declared_dependency_is_documented():
    """docs/DEPENDENCIES.md must mention every library declared in pyproject.toml and superset/requirements.txt."""
    import re
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    doc = (root / "docs" / "DEPENDENCIES.md").read_text(encoding="utf-8").lower()
    py = (root / "pyproject.toml").read_text(encoding="utf-8")
    names = set(re.findall(r'"([A-Za-z0-9_.-]+)\s*(?:[<>=!~][^"]*)?"', py[py.index("dependencies = ["):py.index("[build-system]")]))
    names |= {re.split(r"[<>=!~ ]", ln.strip())[0] for ln in (root / "superset" / "requirements.txt").read_text(encoding="utf-8").splitlines()
              if ln.strip() and not ln.startswith("#")}
    names = {n for n in names if n and n not in ("genomecurator",)}
    assert {"fastapi", "httpx", "apscheduler", "pytest", "pillow"} <= {n.lower() for n in names}   # the parser found the lists
    missing = sorted(n for n in names if n.lower() not in doc)
    assert not missing, f"not documented in docs/DEPENDENCIES.md: {missing} (run scripts/dependencies.py)"


def test_schema_board_maps_fields_to_columns_and_propagates_problems():
    from genomecurator import schema_check as sc
    m = sc.field_map("ncbi", {"length", "virus.organism_name", "location.geographic_location", "host.organism_name", "virus.lineage[].tax_id", "virus.lineage[].name"})
    assert m["length"] == ["genome_length"] and m["virus.organism_name"] == ["genome_name"]
    assert m["location.geographic_location"] == ["geographic_location", "isolation_country"]       # one field feeding two columns
    assert m["virus.lineage[].tax_id"] == ["taxon_lineage_ids"] and m["virus.lineage[].name"] == ["taxon_lineage_names"]
    mp = sc.field_map("pathoplexus", {"hostNameScientific", "hostNameCommon", "accessionVersion"})
    assert mp["hostNameCommon"] == ["host_name"] and mp["accessionVersion"] == ["pathoplexus_accession"]   # fallbacks and renames are found too
    rows = [{"accession": "A.1", "length": 5, "host": {"organism_name": "H"}} for _ in range(4)]
    assert sc.profile([{"virus": {"lineage": [{"tax_id": 1}, {"tax_id": 2}, {"tax_id": 3}]}}])["virus.lineage[].tax_id"]["fill"] == 1.0   # not 300%
    read = {"accession", "length", "host.organism_name"}
    fm = {"length": ["genome_length"], "host.organism_name": ["host_name"], "accession": ["accession"]}
    live_bad = sc.profile([{"accession": "A.1", "seqLength": 5} for _ in range(4)])
    findings = [{"source": "ncbi", "severity": "fail", "kind": "removed", "path": "host.organism_name", "detail": "gone", "hint": ""},
                {"source": "ncbi", "severity": "fail", "kind": "removed", "path": "length", "detail": "gone", "hint": ""}]
    cols = ["genome_length", "host_name", "accession", "clade"]
    b = sc.build_board({"ncbi": {"live": live_bad, "read": read, "fmap": fm}}, findings, set(), cols)
    st = {c["name"]: c["status"] for c in b["local"]["columns"]}
    assert st == {"genome_length": "fail", "host_name": "fail", "accession": "ok", "clade": "internal"}   # every feeder broken -> fail; app-filled -> internal
    assert b["sources"]["ncbi"]["status"] == "fail"
    fl = {f["path"]: f for f in b["sources"]["ncbi"]["fields"]}
    assert fl["host.organism_name"]["status"] == "fail" and fl["host.organism_name"]["group"] == "host" and fl["host.organism_name"]["name"] == "organism_name"
    # a second source that still fills the column makes the problem partial (a warning), not a failure
    b2 = sc.build_board({"ncbi": {"live": live_bad, "read": read, "fmap": fm},
                         "bvbrc": {"live": sc.profile([{"genome_length": 5}]), "read": {"genome_length"}, "fmap": {"genome_length": ["genome_length"]}}}, findings, set(), cols)
    assert {c["name"]: c["status"] for c in b2["local"]["columns"]}["genome_length"] == "warn"
    b3 = sc.build_board({"ncbi": {"live": live_bad, "read": read, "fmap": fm}}, [], {"genome_length"}, cols)
    assert {c["name"]: c["status"] for c in b3["local"]["columns"]}["genome_length"] == "fail"     # missing in the real database


def test_schema_check_does_not_cry_wolf_about_sparse_fields_or_legacy_baselines():
    from genomecurator import schema_check as sc
    base = [{"accession": "A.1", "length": 5, "isolate": {"collection_date": "2020"}} for _ in range(100)]
    read = {"accession", "length", "isolate.collection_date"}
    sh = {"top": ["reports"]}
    # a field that was filled in 1% of records: 200 records without it prove nothing (~13% by pure chance)
    prev = sc.profile([{**r, "biosample": "SAMN1"} if i == 0 else r for i, r in enumerate(base)])
    now = sc.profile(base * 2)
    assert prev["biosample"]["fill"] == 0.01
    assert sc._compare("ncbi", now, read | {"biosample"}, prev, sh, ["reports"], n=200) == []
    # a field that was filled in 94% of records and is gone from 200 is real
    prev94 = sc.profile([{**r, "host": {"organism_name": "H"}} if i < 94 else r for i, r in enumerate(base)])
    f = sc._compare("ncbi", now, read | {"host.organism_name"}, prev94, sh, ["reports"], n=200)
    assert any(x["kind"] == "removed" and x["path"] == "host.organism_name" for x in f)
    # a baseline written by older code with an impossible fill (500%) must not trigger a bogus "dropped" warning
    legacy = {p: dict(v) for p, v in sc.profile(base).items()}
    legacy["length"]["fill"] = 5.0
    assert sc._compare("ncbi", sc.profile(base), read, legacy, sh, ["reports"], n=100) == []


def test_old_schema_reports_still_render_after_the_code_grows():
    """A board saved before the simple view existed has no column grouping; reading it must not break the Schema tab."""
    from genomecurator import schema_check as sc
    old = {"board": {"local": {"columns": [{"name": "genome_length", "status": "ok", "feeders": []}, {"name": "weird_new_col", "status": "internal"}]}}}
    cols = sc.upgrade_report(old)["board"]["local"]["columns"]
    assert cols[0]["group"] == "유전체" and cols[0]["group_order"] == 2
    assert cols[1]["group"] == "기타" and cols[1]["feeders"] == [] and cols[1]["msg"] == ""
    assert sc.upgrade_report({"verdict": "ok"}) == {"verdict": "ok"}        # a report without a board is untouched


def test_completeness_classes_and_sql_view_agree(tmp_path, monkeypatch):
    from genomecurator import completeness as C, config, db
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    db._recorded.discard(str(config.DB_PATH))
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'V','V',9),(2,'NoRef','N',8)")
    rows = [("REF", 1000, "Complete", "NC_000001"), ("c-short", 500, "Complete", None),          # NCBI's Complete is trusted even when short
            ("edge90", 900, "Partial", None), ("below", 899, "Partial", None), ("long", 1000, "Partial", None),
            ("nostatus", 950, None, None), ("gone", 950, "Removed", None)]
    for i, (acc, ln, st, rs) in enumerate(rows):
        con.execute("INSERT INTO genome(genome_id,watch_id,accession,genome_length,genome_status,refseq_accessions) VALUES(?,?,?,?,?,?)", (str(i), 1, acc, ln, st, rs))
    con.execute("INSERT INTO genome(genome_id,watch_id,accession,genome_length,genome_status) VALUES('x',2,'noref',5000,'Partial')")   # no RefSeq: cannot be Near complete
    con.commit()
    ref = C.reference_map(con)
    assert ref == {1: ("REF", 1000)}
    for pct in (90.0, 95.0):                       # the view reads the threshold stored in app_meta, so both paths must follow the setting
        monkeypatch.setattr(config, "NEAR_COMPLETE_PCT", pct)
        db._recorded.discard(str(config.DB_PATH))
        con = connect()
        sql = {r[0]: r[1] for r in con.execute("SELECT accession, completeness FROM v_genome_analysis")}
        for acc, ln, st, wid in con.execute("SELECT accession, genome_length, genome_status, watch_id FROM genome"):
            assert sql[acc] == C.classify(st, ln, ref.get(wid, (None, None))[1], pct)[0], (acc, pct)
        con.close()
    assert sql["edge90"] == "Partial" and C.classify("Partial", 900, 1000, 90)[0] == "Near complete"   # 95% setting is stricter than 90%
    assert sql["c-short"] == "Complete" and sql["noref"] == "Partial" and sql["gone"] == "Removed"
    assert "판정하지 않습니다" in C.describe(None) and "NC_000001" in C.describe(("NC_000001", 1000)) and "90%" in C.describe(("NC_000001", 1000), 90)


def test_reconcile_host_ignores_a_parenthesised_common_name():
    from genomecurator.reconcile import classify
    assert classify("host", "Sus scrofa", "Sus scrofa (pig)") == "cosmetic"
    assert classify("host", "Pteropus hypomelanus", "Pteropus hypomelanus (Island flying-fox)") == "cosmetic"
    assert classify("host", "Homo sapiens", "Homo sapiens; 32; male") == "cosmetic"
    assert classify("host", "Pteropus vampyrus", "Pteropus giganteus (bat)") == "conflict"


def test_resolve_suggestions_follow_the_stated_rules():
    from genomecurator import resolve as R
    g = lambda *pairs: [{"value": v, "sources": s} for v, s in pairs]
    # same fact, different writing: the most specific date, ISO preferred at equal precision
    s = R.suggest("collection_date", g(("2022-07", ["ncbi"]), ("2022-07-11", ["bvbrc", "pathoplexus"])))
    assert s["value"] == "2022-07-11" and s["confidence"] == "high"
    s = R.suggest("collection_date", g(("2022-06-19", ["ncbi"]), ("19-Jun-2022", ["bvbrc"])))
    assert s["value"] == "2022-06-19" and s["confidence"] == "high"
    # an impossible date (in the future, or before 1900) is never suggested, even when it is the priority source
    s = R.suggest("collection_date", g(("2099-05-30", ["ncbi"]), ("2019-05-30", ["bvbrc"])))
    assert s["value"] == "2019-05-30" and "불가능" in s["reason"]
    s = R.suggest("collection_date", g(("1822-05-30", ["pathoplexus"]), ("2022-05-30", ["bvbrc"])))
    assert s["value"] == "2022-05-30" and "불가능" in s["reason"]
    # ... but 1922 is a *possible* date (historical specimens exist), so two sources that disagree on it leave the decision to the person
    s = R.suggest("collection_date", g(("1922-05-30", ["ncbi"]), ("2022-05-30", ["bvbrc"])))
    assert s["confidence"] == "low"
    # a majority wins and is MEDIUM; a minority date that disagrees is outvoted
    s = R.suggest("collection_date", g(("2022-07-25", ["ncbi", "pathoplexus"]), ("2023", ["bvbrc"])))
    assert s["value"] == "2022-07-25" and s["confidence"] == "medium" and "3개 소스 중 2개" in s["reason"]
    # length / version: NCBI, with the reason (the sequence is downloaded from NCBI and verified against it)
    s = R.suggest("genome_length", g(("197179", ["ncbi", "pathoplexus"]), ("197157", ["bvbrc"])), {"pathoplexus": "same", "bvbrc": "diff"})
    assert s["value"] == "197179" and s["confidence"] == "high" and "Pathoplexus" in s["reason"]
    # host: the scientific name is the canonical writing
    s = R.suggest("host_name", g(("Sus scrofa", ["ncbi"]), ("Sus scrofa (pig)", ["bvbrc"])))
    assert s["value"] == "Sus scrofa" and s["confidence"] == "high"
    # two sources, really different values: no basis to choose -> LOW and the person is told to look
    s = R.suggest("isolation_source", g(("oronasopharynx", ["ncbi"]), ("throat swab", ["bvbrc"])))
    assert s["confidence"] == "low" and "직접 확인" in s["reason"] and s["value"] == "oronasopharynx"
    # more specific wording for a pure rewording
    assert R.suggest("isolation_source", g(("urine", ["ncbi"]), ("urine sample", ["bvbrc"])))["value"] == "urine sample"


def test_resolve_contained_wording_is_medium_not_high():
    from genomecurator import resolve as R
    g = lambda *pairs: [{"value": v, "sources": s} for v, s in pairs]
    assert R.suggest("isolation_source", g(("liver, spleen", ["ncbi"]), ("Liver/Spleen", ["bvbrc"])))["confidence"] == "high"      # only punctuation / case differ
    s = R.suggest("isolation_source", g(("urine", ["ncbi"]), ("urine of Pteropus vampyrus", ["bvbrc"])))
    assert s["confidence"] == "medium" and s["value"] == "urine of Pteropus vampyrus" and "포함" in s["reason"]


def _resolve_fixture(tmp_path, monkeypatch):
    import json
    from genomecurator import config, db
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    db._recorded.discard(str(config.DB_PATH))
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'V','V',9)")
    src = {"ncbi": {"A1": {"collection_date": "2022-07", "host_name": "Chiroptera", "genome_length": 100}},
           "bvbrc": {"A1": {"collection_date": "2022-07-11", "host_name": "bat", "genome_length": 90}},
           "pathoplexus": {"A1": {"collection_date": "2022-07-11", "host_name": "Chiroptera"}}}
    for s, recs in src.items():
        for a, r in recs.items():
            con.execute("INSERT INTO source_record(source,accession,watch_id,payload_json) VALUES(?,?,1,?)", (s, a, json.dumps(r)))
    merged = {"collection_date": "2022-07-11", "host_name": "Chiroptera", "genome_length": 100}      # what the priority merge gave
    con.execute("INSERT INTO changeset(id,watch_id,status) VALUES(1,1,'READY')")
    con.execute("INSERT INTO changeset_item(id,changeset_id,accession,kind,diff_json,new_record_json,provenance_json,fingerprint,decision) VALUES(7,1,'A1','new','{}',?,?,'x','pending')",
                (json.dumps(merged), json.dumps({"collection_date": "pathoplexus", "host_name": "pathoplexus", "genome_length": "ncbi"})))
    con.commit()
    return con, src


def test_resolve_save_applies_to_the_item_and_survives_the_next_merge(tmp_path, monkeypatch):
    import json
    from genomecurator import curate, resolve as R
    con, src = _resolve_fixture(tmp_path, monkeypatch)
    view = R.item_view(con, 1, 7)
    fields = {f["field"]: f for f in view["fields"]}
    assert set(fields) == {"collection_date", "host_name", "genome_length"} and fields["genome_length"]["suggestion"]["value"] == "100"
    assert fields["collection_date"]["suggestion"]["value"] == "2022-07-11" and fields["collection_date"]["suggestion"]["confidence"] == "high"
    assert fields["host_name"]["suggestion"]["value"] == "Chiroptera" and fields["host_name"]["suggestion"]["confidence"] == "medium"    # 2 of 3 sources
    # the person picks BV-BRC's host and types their own date
    con.execute("BEGIN")
    out = R.save(con, 1, "A1", [{"field": "host_name", "source": "bvbrc", "value": "bat"}, {"field": "collection_date", "source": "custom", "value": "2022-07-12"}])
    con.commit()
    assert out == {"saved": 2, "reset": 0, "items": 1}
    row = con.execute("SELECT * FROM changeset_item WHERE id=7").fetchone()
    rec, prov = json.loads(row["new_record_json"]), json.loads(row["provenance_json"])
    assert rec["host_name"] == "bat" and prov["host_name"] == "bvbrc" and rec["collection_date"] == "2022-07-12" and prov["collection_date"] == "manual"
    assert row["fingerprint"] != "x" and json.loads(row["diff_json"])["host_name"] == [None, "bat"]
    # a later check merges again from scratch: the choices must still win
    merged = curate.merge({"ncbi": {"A1": src["ncbi"]["A1"]}, "bvbrc": {"A1": src["bvbrc"]["A1"]}, "pathoplexus": {"A1": src["pathoplexus"]["A1"]}}, R.load_choices(con, 1))
    assert merged["A1"][0]["host_name"] == "bat" and merged["A1"][1]["host_name"] == "bvbrc" and merged["A1"][0]["collection_date"] == "2022-07-12"
    # the chosen source's CURRENT value is used (it changed), and if that source loses the field the priority merge stands
    changed = {"ncbi": {"A1": src["ncbi"]["A1"]}, "bvbrc": {"A1": {**src["bvbrc"]["A1"], "host_name": "Bat"}}, "pathoplexus": {"A1": src["pathoplexus"]["A1"]}}
    assert curate.merge(changed, R.load_choices(con, 1))["A1"][0]["host_name"] == "Bat"
    gone = {"ncbi": {"A1": src["ncbi"]["A1"]}, "bvbrc": {"A1": {"genome_length": 90}}, "pathoplexus": {"A1": src["pathoplexus"]["A1"]}}
    assert curate.merge(gone, R.load_choices(con, 1))["A1"][0]["host_name"] == "Chiroptera"
    # decided fields no longer count as open disagreements
    con.executemany("INSERT INTO source_conflict(watch_id,accession,field,kind,src_a,val_a,src_b,val_b) VALUES(1,'A1',?,?,'ncbi','a','bvbrc','b')",
                    [("host_name", "conflict"), ("collection_date", "cosmetic"), ("genome_length", "conflict")])
    con.commit()
    ag = curate.source_agreement(1)["A1"]
    assert ag["conflict"] == {"genome_length"} and ag["cosmetic"] == set() and ag["chosen"] == {"host_name", "collection_date"}
    # back to the default for one field
    con.execute("BEGIN")
    R.save(con, 1, "A1", [{"field": "host_name", "source": None}])
    con.commit()
    assert json.loads(con.execute("SELECT new_record_json FROM changeset_item WHERE id=7").fetchone()[0])["host_name"] == "Chiroptera"
    assert set(R.load_choices(con, 1)["A1"]) == {"collection_date"}
    import pytest
    for bad in ({"field": "nope", "source": "ncbi"}, {"field": "host_name", "source": "custom", "value": " "}, {"field": "genome_length", "source": "pathoplexus"}):
        with pytest.raises(ValueError):
            R.save(con, 1, "A1", [bad])


def test_resolve_auto_pick_applies_medium_and_leaves_low(tmp_path, monkeypatch):
    import json
    from genomecurator import resolve as R
    con, _ = _resolve_fixture(tmp_path, monkeypatch)
    con.execute("DELETE FROM source_record WHERE source='pathoplexus'")          # two sources only: host is now a 1-vs-1 disagreement = LOW
    con.commit()
    pv = R.preview_auto(con, 1, "medium")
    assert pv["counts"]["low"] == 1 and pv["counts"]["high"] == 2 and pv["fields"] == 2         # dates and length are certain, the host is left to the person
    res = R.apply_auto(con, 1, "medium")
    con.commit()
    assert res["applied"] == 2 and set(R.load_choices(con, 1)["A1"]) == {"collection_date", "genome_length"}
    assert all(r[0] == 1 for r in con.execute("SELECT auto FROM field_choice"))             # marked as automatic
    assert R.preview_auto(con, 1, "medium")["fields"] == 0                                    # nothing left to pick automatically
    assert "field_choice" in __import__("genomecurator.backup", fromlist=["x"]).KEEP_TABLES   # decisions are part of backups


def _jobs_env(tmp_path, monkeypatch):
    from genomecurator import config, db, jobs
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    db.connect().close()
    jobs.live.clear()
    return jobs


def test_job_lifecycle_is_recorded_with_log_error_and_result(tmp_path, monkeypatch):
    jobs = _jobs_env(tmp_path, monkeypatch)

    def ok(p):
        p(40, "halfway")
        return {"redirect": "/x", "result": {"found": 3}}

    def boom(p):
        p(10, "starting")
        raise RuntimeError("source down")
    a = jobs.run("aaaa0001", "good job", "check", ok)
    b = jobs.run("bbbb0002", "bad job", "apply", boom)
    assert a["status"] == "done" and a["pct"] == 100 and a["redirect"] == "/x" and a["result"] == {"found": 3}
    assert b["status"] == "error" and "source down" in b["error"] and any(x["level"] == "error" for x in b["log"])
    jobs.live.clear()                                        # as after the live entry expired (or a reload): the table is the record
    assert jobs.get("aaaa0001")["result"] == {"found": 3} and any("halfway" in x["msg"] for x in jobs.get("aaaa0001")["log"])
    assert [r["id"] for r in jobs.list_jobs(status="error")] == ["bbbb0002"]
    assert jobs.get("nope") is None


def test_running_job_is_interrupted_by_a_restart_and_summary_counts(tmp_path, monkeypatch):
    jobs = _jobs_env(tmp_path, monkeypatch)
    import threading, time
    gate, started = threading.Event(), threading.Event()

    def slow(p):
        p(20, "working")
        started.set()
        gate.wait(5)
    jid = jobs.start(slow, "slow one", "fasta")
    assert started.wait(5)
    s = jobs.summary()
    assert [r["id"] for r in s["running"]] == [jid]
    gate.set()
    for _ in range(100):
        if jobs.live[jid]["status"] != "running":
            break
        time.sleep(0.02)
    # a row the previous process left as running
    from genomecurator.db import connect
    con = connect()
    con.execute("INSERT INTO job(id,kind,label,status,pct,msg,started) VALUES('dead0001','check','died','running',30,'x','2026-01-01T00:00:00')")
    con.commit()
    con.close()
    assert jobs.mark_interrupted() == 1
    d = jobs.get("dead0001")
    assert d["status"] == "interrupted" and "다시 시작" in d["error"]
    assert jobs.summary()["failed"] >= 1


def test_jobs_api_and_server_side_add_and_lookup(tmp_path, monkeypatch):
    jobs = _jobs_env(tmp_path, monkeypatch)
    import time
    from fastapi.testclient import TestClient
    from genomecurator import curate, lookup, web
    calls = []
    monkeypatch.setattr(curate, "add_watch", lambda name, taxon=None, filters=None: calls.append((name, taxon)) or 1)
    monkeypatch.setattr(lookup, "lookup_virus", lambda q: {"sources": {}, "add": {"taxon": "T", "taxon_id": 1}, "suggestions": [], "note": q})
    c = TestClient(web.create_app(scheduler=False))

    def wait(jid):
        for _ in range(100):
            j = c.get(f"/api/job/{jid}").json()
            if j["done"]:
                return j
            time.sleep(0.05)
        raise AssertionError("job did not finish")
    j = wait(c.post("/watch", data={"name": "Mpox", "taxon": "Monkeypox virus"}).json()["job"])
    assert j["status"] == "done" and calls == [("Mpox", "Monkeypox virus")] and j["redirect"].startswith("/?msg=")
    j = wait(c.post("/api/lookup/start?q=nipah").json()["job"])
    assert j["result"]["note"] == "nipah" and j["kind"] == "lookup"
    assert c.post("/api/lookup/start?q=").status_code == 400
    monkeypatch.setattr(curate, "add_watch", lambda *a, **k: (_ for _ in ()).throw(ValueError("no such taxon")))
    j = wait(c.post("/watch", data={"name": "zzz"}).json()["job"])
    assert j["status"] == "error" and "no such taxon" in j["error"]       # the failure is visible, not swallowed
    assert {r["kind"] for r in c.get("/api/jobs").json()} == {"watch_add", "lookup"}
    assert c.get("/api/jobs/summary").json()["failed"] >= 1
    assert c.get("/jobs").status_code == 200 and c.get("/history?q=job%3D" + j["id"]).status_code == 200


def test_analysis_page_and_superset_start_reports_not_installed(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from genomecurator import analysis, config, db, web
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    monkeypatch.setattr(analysis, "HOME", tmp_path / "nohome")           # not initialised
    monkeypatch.setattr(analysis, "listening", lambda: False)
    db.connect().close()
    c = TestClient(web.create_app(scheduler=False))
    r = c.get("/analysis")
    assert r.status_code == 200 and "아직 설치되지 않았습니다" in r.text
    assert c.get("/api/superset/status").json()["installed"] is False
    import time
    jid = c.post("/api/superset/start").json()["job"]
    for _ in range(100):
        j = c.get(f"/api/job/{jid}").json()
        if j["done"]:
            break
        time.sleep(0.05)
    assert j["status"] == "error" and "설치" in j["error"]              # the reason is shown, nothing is launched
    monkeypatch.setattr(analysis, "listening", lambda: True)
    assert analysis.start() == f"http://127.0.0.1:{analysis.port()}/"    # already running: nothing to start


def test_resolve_auto_for_selected_items_only(tmp_path, monkeypatch):
    from genomecurator import resolve as R
    con, _ = _resolve_fixture(tmp_path, monkeypatch)
    item = con.execute("SELECT id FROM changeset_item WHERE accession='A1'").fetchone()[0]
    all_ = R.preview_auto(con, 1, "low")
    assert all_["accessions"] == 1
    assert R.preview_auto(con, 1, "low", [item + 1000])["items"] == 0          # a row that is not selected / does not exist: nothing
    assert R.preview_auto(con, 1, "low", [])["items"] == 0                      # nothing selected: nothing (never "everything")
    assert R.preview_auto(con, 1, "low", [item])["fields"] == all_["fields"]
    res = R.apply_auto(con, 1, "low", [item])
    con.commit()
    assert res["applied"] == all_["fields"] and set(R.load_choices(con, 1)["A1"])


def test_resolve_auto_api_takes_selected_items(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from genomecurator import web
    con, _ = _resolve_fixture(tmp_path, monkeypatch)
    item = con.execute("SELECT id FROM changeset_item WHERE accession='A1'").fetchone()[0]
    con.close()
    c = TestClient(web.create_app(scheduler=False))
    assert c.post("/api/changeset/1/resolve_auto/preview", json={"items": [item + 1000], "minimum": "low"}).json()["items"] == 0
    pv = c.post("/api/changeset/1/resolve_auto/preview", json={"items": [item], "minimum": "low"}).json()
    assert pv["items"] == 1 and pv["fields"] >= 1
    assert c.post("/changeset/1/resolve_auto", json={"items": [item + 1000], "minimum": "low"}).json()["applied"] == 0
    assert c.post("/changeset/1/resolve_auto", json={"items": [item], "minimum": "low"}).json()["applied"] == pv["fields"]
