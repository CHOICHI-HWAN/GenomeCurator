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


def test_schema_traffic_light_per_source(tmp_path, monkeypatch):
    import json
    from genomecurator import config, schema_check
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    t = schema_check.traffic()
    assert {s["status"] for s in t["sources"].values()} == {"unknown"} and list(t["sources"]) == ["ncbi", "pathoplexus", "bvbrc"]
    rep = {"mode": "daily", "checked_at": "2026-10-19T03:00:00", "verdict": "update", "findings": [
        {"source": "bvbrc", "severity": "warn", "path": "a", "detail": "x"}, {"source": "pathoplexus", "severity": "fail", "path": "b", "detail": "gone"},
        {"source": "pathoplexus", "severity": "warn", "path": "c", "detail": "y"}, {"source": "local", "severity": "fail", "path": "-", "detail": "z"}]}
    con = connect()
    con.execute("INSERT INTO schema_report(checked_at, report_json) VALUES(?,?)", (rep["checked_at"], json.dumps(rep)))
    con.commit()
    con.close()
    t = schema_check.traffic()
    assert {k: v["status"] for k, v in t["sources"].items()} == {"ncbi": "ok", "pathoplexus": "fail", "bvbrc": "warn"}
    assert t["sources"]["pathoplexus"]["fail"] == 1 and "gone" in t["sources"]["pathoplexus"]["msg"]


def _qc_env(tmp_path, monkeypatch, n_N=0):
    """A watch with two stored genomes (A1, X2) and a fake `nextclade` on PATH: it calls A1 clade Ib and X2 clade IIb."""
    import os
    import stat
    from genomecurator import config, db
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    nc = bindir / ("nextclade.py" if os.name == "nt" else "nextclade")
    nc.write_text("""#!/usr/bin/env python3
import sys, json, os
a = sys.argv[1:]
if a[0] == "--version":
    print("nextclade 9.9.9")
elif a[:2] == ["dataset", "list"]:
    print(json.dumps([{"version": {"tag": "2030-01-01--00-00-00Z"}}]))
elif a[:2] == ["dataset", "get"]:
    d = a[a.index("--output-dir") + 1]
    os.makedirs(d)
    open(os.path.join(d, "pathogen.json"), "w").write(json.dumps({"version": {"tag": "2030-01-01--00-00-00Z"}}))
elif a[0] == "run":
    out = a[a.index("--output-tsv") + 1]
    names = [l[1:].strip() for l in open(a[-1]) if l.startswith(">")]
    with open(out, "w") as f:
        f.write("index\\tseqName\\tclade\\tlineage\\tqc.overallScore\\tqc.overallStatus\\tcoverage\\ttotalMissing\\tqc.privateMutations.status\\terrors\\n")
        for i, n in enumerate(names):
            f.write(f"{i}\\t{n}\\t{'IIb' if n.startswith('X') else 'Ib'}\\tA.4\\t{5 * i}\\t{'bad' if n.startswith('X') else 'good'}\\t0.99\\t3\\t{'mediocre' if n.startswith('X') else 'good'}\\t\\n")
""")
    nc.chmod(nc.stat().st_mode | stat.S_IEXEC)
    if os.name == "nt":      # Windows runs no #! scripts: a .cmd wrapper that shutil.which finds through PATHEXT
        import sys
        (bindir / "nextclade.cmd").write_text(f'@"{sys.executable}" "{nc}" %*\r\n')
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    con = db.connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'Mpox','Orthopoxvirus monkeypox',3431483)")
    for i, (acc, seq) in enumerate([("A1", "ACGT" * 50), ("X2", "ACGT" * 40 + "N" * 40)]):
        (tmp_path / "genomes").mkdir(exist_ok=True)
        f = tmp_path / "genomes" / f"{acc}.fna"
        f.write_text(f">{acc} test\n{seq}\n")
        con.execute("INSERT INTO genome(genome_id,watch_id,accession,genome_length,clade) VALUES(?,?,?,?,?)", (f"g{i}", 1, acc, len(seq), "Ib" if acc == "A1" else "IIa"))
        con.execute("INSERT INTO genome_sequence(accession,genome_id,md5,file_path) VALUES(?,?,?,?)", (acc, f"g{i}", f"md5{acc}", f"genomes/{acc}.fna"))
    con.commit()
    con.close()


def test_qc_seqstats_verdicts(tmp_path):
    from genomecurator import qc
    f = tmp_path / "a.fna"
    f.write_text(">a\n" + "ACGT" * 100 + "\n")
    assert qc.seqstats(f, 400)["status"] == "good"
    assert qc.seqstats(f, 399)["status"] == "bad"                       # length differs from the record
    f.write_text(">a\n" + "ACGT" * 95 + "N" * 20 + "\n")                  # 5% N
    assert qc.seqstats(f)["status"] == "bad"
    f.write_text(">a\n" + "ACGT" * 99 + "NNNN\n")                          # 1% N
    assert qc.seqstats(f)["status"] == "mediocre"
    f.write_text(">a\n" + "ACGT" * 100 + "R" + "\n")                       # a single ambiguity code in 401: 0.25%
    assert qc.seqstats(f)["status"] == "mediocre"
    f.write_text(">a\n")
    assert qc.seqstats(f)["status"] == "bad"


def test_qc_run_stores_results_and_skips_unchanged(tmp_path, monkeypatch):
    from genomecurator import qc
    _qc_env(tmp_path, monkeypatch)
    assert qc.dataset_for("Orthopoxvirus monkeypox") == "nextstrain/mpox/all-clades" and qc.dataset_for("Nipah virus") is None
    st = qc.dataset_status()[0]
    assert st["state"] == "missing" and st["latest"] == "2030-01-01--00-00-00Z"          # nothing downloaded yet
    r0 = qc.run(1)
    assert r0["nextclade"] == 0 and r0["stats"] == 2 and any("내려받아지지" in n for n in r0["notes"])       # stats only, and it says why
    up = qc.dataset_update("nextstrain/mpox/all-clades")
    assert up["after"] == "2030-01-01--00-00-00Z" and qc.dataset_status()[0]["state"] == "current"
    r = qc.run(1)
    assert r["nextclade"] == 2 and r["skipped"] == 2                                       # stats were already there
    pg = qc.per_genome(1)
    assert pg["A1"]["qc"] == "good" and pg["A1"]["nc_clade"] == "Ib"
    assert pg["X2"]["qc"] == "bad" and pg["X2"]["nc_clade"] == "IIb" and pg["X2"]["stats_status"] == "bad"
    again = qc.run(1)
    assert again["nextclade"] == 0 and again["stats"] == 0                                 # nothing changed: nothing recomputed
    assert qc.run(1, force=True)["nextclade"] == 2
    s = qc.summary(1)
    assert s["clade_diff"] == 1 and s["clade_examples"][0]["accession"] == "X2"            # source says IIa, Nextclade says IIb
    assert qc.nextclade_source(1) == {"A1": {"clade": "Ib", "lineage": "A.4"}, "X2": {"clade": "IIb", "lineage": "A.4"}}


def test_nextclade_enters_curation_as_a_source(tmp_path, monkeypatch):
    from genomecurator import config, curate, qc, reconcile, resolve
    _qc_env(tmp_path, monkeypatch)
    qc.dataset_update("nextstrain/mpox/all-clades")
    qc.run(1)
    ps = {"ncbi": {"A1": {"genome_length": 200, "clade": "Ib"}, "X2": {"genome_length": 200, "clade": "IIa"}, "Y3": {"genome_length": 5}},
          "bvbrc": {}, "pathoplexus": {}, "nextclade": qc.nextclade_source(1)}
    rep, rows = reconcile.analyze(ps)
    assert rep["coverage"]["sources"] == {"ncbi": 3}                                       # Nextclade is not a database: no coverage numbers
    assert [(r[0], r[1], r[2], r[3], r[5]) for r in rows] == [("X2", "clade", "conflict", "ncbi", "nextclade")]
    m = curate.merge(ps)
    assert m["A1"][0]["clade"] == "Ib" and m["X2"][0]["clade"] == "IIa" and m["X2"][1]["clade"] == "ncbi"       # a database wins by default
    ps["ncbi"]["Y3"].pop("genome_length")
    ps["ncbi"]["X2"].pop("clade")
    assert curate.merge(ps)["X2"][0]["clade"] == "IIb" and curate.merge(ps)["X2"][1]["clade"] == "nextclade"   # but Nextclade fills a blank
    fs = resolve.analyse({"ncbi": {"clade": "IIa"}, "nextclade": {"clade": "IIb"}}, {"clade": "IIa"}, {"clade": "ncbi"}, {})
    assert fs[0]["suggestion"]["value"] == "IIb" and fs[0]["suggestion"]["confidence"] == "medium"
    fs = resolve.analyse({"ncbi": {"clade": "IIb"}, "bvbrc": {"clade": "IIa"}, "nextclade": {"clade": "IIb"}}, {}, {}, {})
    assert fs[0]["suggestion"]["value"] == "IIb" and fs[0]["suggestion"]["confidence"] == "high"


def test_qc_api_and_pages(tmp_path, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from genomecurator import web
    _qc_env(tmp_path, monkeypatch)
    c = TestClient(web.create_app(scheduler=False))
    assert c.get("/qc").status_code == 200
    d = c.get("/api/qc/datasets").json()
    assert d["tool"] == "9.9.9" and d["rows"][0]["state"] == "missing"
    assert c.post("/api/qc/dataset?path=evil/path").status_code == 400                   # only datasets of registered viruses

    def wait(jid):
        for _ in range(200):
            j = c.get(f"/api/job/{jid}").json()
            if j["done"]:
                return j
            time.sleep(0.05)
        raise AssertionError("job did not finish")
    assert wait(c.post("/api/qc/dataset?path=nextstrain/mpox/all-clades").json()["job"])["status"] == "done"
    j = wait(c.post("/api/qc/run?watch=1").json()["job"])
    assert j["status"] == "done" and j["result"]["nextclade"] == 2 and j["kind"] == "qc"
    assert "QC" in c.get("/qc").text and c.get("/genomes?watch=1").status_code == 200
    assert c.post("/api/qc/run?watch=99").status_code == 404


def test_export_helpers():
    from genomecurator import export
    assert export.std_date("2023") == "2023-XX-XX" and export.std_date("2023-06") == "2023-06-XX" and export.std_date("2023-06-15") == "2023-06-15"
    assert export.std_date("Mar-2026") == "2026-03-XX" and export.std_date("") == "" and export.std_date("2023-01/2023-03") == "" and export.std_date("garbage") == ""
    assert [export.region_of(c) for c in ("USA", "Democratic Republic of the Congo", "Germany", "Viet Nam", "Brazil", "New Zealand", "Atlantis", "")] == \
        ["North America", "Africa", "Europe", "Asia", "South America", "Oceania", "", ""]
    assert export._geo({"geographic_location": "USA, NY, New York City", "isolation_country": "USA"}) == ("NY", "New York City")
    assert export._geo({"geographic_location": "Uganda: Kampala", "isolation_country": "Uganda"}) == ("Kampala", "")
    assert export._geo({"geographic_location": "Germany", "isolation_country": "Germany"}) == ("", "")


def test_export_writes_aligned_metadata_and_fasta(tmp_path, monkeypatch):
    import hashlib
    import json
    from genomecurator import export, qc
    _qc_env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    con = connect()
    con.execute("UPDATE genome SET isolation_country='USA', geographic_location='USA, NY, New York City', collection_date='2023-06', host_name='Homo sapiens', "
                "submitter='Lab X', accession_version=1, pathoplexus_accession='PP_000TEST.2', release_date='2024-01-02T03:04:05Z' WHERE accession='A1'")
    con.execute("UPDATE genome SET pathoplexus_accession='PP_000TX2.1' WHERE accession='X2'")
    con.commit()
    con.close()
    qc.dataset_update("nextstrain/mpox/all-clades")
    qc.run(1)
    out = tmp_path / "out"
    r = export.write(1, out, "nextstrain")
    assert r["rows"] == 2 and r["sequences"] == 2
    lines = (out / "metadata.tsv").read_text(encoding="utf-8").splitlines()
    cols = lines[0].split("\t")
    assert cols == export.NEXTSTRAIN_COLUMNS
    rows = [dict(zip(cols, ln.split("\t"))) for ln in lines[1:]]
    ids = [x["PPX_accession"] for x in rows]
    fa_ids = [ln[1:] for ln in (out / "sequences.fasta").read_text(encoding="utf-8").splitlines() if ln.startswith(">")]
    assert ids == fa_ids == ["PP_000TEST", "PP_000TX2"]                      # same ids, same order; the id is the PPX accession the workflow keys on
    assert rows[0]["INSDC_accession_version"] == "A1.1" and rows[1]["INSDC_accession_version"] == "X2"
    a = rows[0]
    assert (a["date"], a["country"], a["region"], a["division"], a["location"], a["host"], a["date_released"]) == \
        ("2023-06-XX", "USA", "North America", "NY", "New York City", "Homo sapiens", "2024-01-02")
    assert (a["PPX_accession"], a["PPX_accession_version"], a["dataUseTerms"], a["institution"], a["clade"], a["QC_rare_mutations"], a["missing_data"]) == \
        ("PP_000TEST", "PP_000TEST.2", "OPEN", "Lab X", "Ib", "good", "3")
    assert rows[1]["clade"] == "IIa" and rows[1]["QC_rare_mutations"] == "mediocre"          # the curated clade wins over Nextclade's IIb
    mf = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert mf["files"]["metadata.tsv"]["sha256"] == hashlib.sha256((out / "metadata.tsv").read_bytes()).hexdigest()
    assert mf["files"]["sequences.fasta"]["sha256"] == hashlib.sha256((out / "sequences.fasta").read_bytes()).hexdigest()
    cfg = (out / "nextstrain-config.yaml").read_text(encoding="utf-8")
    assert mf["nextclade"]["datasets"] and "use_remote_files: false" in cfg and "strain_id_field" not in cfg        # the workflow's default id is PPX_accession
    assert "--configfile defaults/mpxv/config.yaml" in cfg and 'metadata: "' in cfg
    seq = "".join(ln for ln in (out / "sequences.fasta").read_text(encoding="utf-8").splitlines()[1:] if not ln.startswith(">"))
    assert seq.startswith("ACGT" * 10)
    # filters
    assert [x["PPX_accession"] for x in export.table(1, "nextstrain", qc_not=["bad"])[1]] == ["PP_000TEST"]             # X2's Nextclade QC is bad
    assert export.table(1, "nextstrain", min_length=201)[1] == []
    assert len(export.table(1, "plain")[1]) == 2 and "nextclade_clade" in export.table(1, "plain")[0]
    # a genome without a Pathoplexus accession cannot be identified for the workflow: it is left out and the preview says so
    con = connect()
    con.execute("UPDATE genome SET pathoplexus_accession=NULL WHERE accession='X2'")
    con.commit()
    con.close()
    pv = export.preview(1, "nextstrain")
    assert pv["selected"] == 1 and pv["left_out_no_id"] == 1 and any("PPX_accession" in w for w in pv["warnings"])
    assert export.preview(1, "plain")["selected"] == 2
    # straight into a nextstrain/mpox checkout: relative paths (they also work inside the Docker runtime)
    ph = tmp_path / "mpox" / "phylogenetic"
    ph.mkdir(parents=True)
    (ph / "Snakefile").write_text("")
    r2 = export.write(1, None, "nextstrain", nextstrain_dir=ph)
    assert (ph / "data" / "genomecurator" / "metadata.tsv").exists() and r2["rows"] == 1
    cfg2 = (ph / "build-configs" / "genomecurator" / "config.yaml").read_text(encoding="utf-8")
    assert 'metadata: "data/genomecurator/metadata.tsv"' in cfg2 and "build-configs/genomecurator/config.yaml" in cfg2
    import pytest
    with pytest.raises(ValueError):
        export.write(1, None, "nextstrain", nextstrain_dir=tmp_path)


def test_export_api(tmp_path, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from genomecurator import qc, web
    _qc_env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    con = connect()
    con.execute("UPDATE genome SET pathoplexus_accession='PP_A1.1' WHERE accession='A1'")
    con.execute("UPDATE genome SET pathoplexus_accession='PP_X2.1' WHERE accession='X2'")
    con.commit()
    con.close()
    qc.dataset_update("nextstrain/mpox/all-clades")
    qc.run(1)
    c = TestClient(web.create_app(scheduler=False))
    assert c.get("/export").status_code == 200 and "nextstrain" in c.get("/api/export/profiles").json()
    pv = c.get("/api/export/1/preview?profile=nextstrain&qc_not=bad").json()
    assert pv["selected"] == 1 and pv["genomes"] == 2
    t = c.get("/api/export/1/metadata.tsv?profile=nextstrain")
    f = c.get("/api/export/1/sequences.fasta?profile=nextstrain")
    assert t.status_code == 200 and f.status_code == 200 and "tab-separated" in t.headers["content-type"]
    assert [ln.split("\t")[0] for ln in t.text.splitlines()[1:]] == [ln[1:] for ln in f.text.splitlines() if ln.startswith(">")] == ["PP_A1", "PP_X2"]
    assert c.get("/api/export/1/metadata.tsv?profile=nope").status_code == 400 and c.get("/api/export/99/preview").status_code == 404
    jid = c.post("/api/export/1/build?profile=nextstrain").json()["job"]
    for _ in range(100):
        j = c.get(f"/api/job/{jid}").json()
        if j["done"]:
            break
        time.sleep(0.05)
    assert j["status"] == "done" and j["result"]["rows"] == 2 and (tmp_path / "exports").exists()
    assert c.get("/api/exports").json()[0]["rows"] == 2


def _priv_files(rows, seqs):
    cols = ["PPX_accession", "strain", "date", "country", "division", "host", "clade"]
    meta = "\t".join(cols) + "\n" + "".join("\t".join(r.get(c, "") for c in cols) + "\n" for r in rows)
    return meta, "".join(f">{k}\n{v}\n" for k, v in seqs.items())


def test_private_check_rejects_structural_problems_and_normalises(tmp_path, monkeypatch):
    from genomecurator import private
    _qc_env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    con = connect()
    con.execute("UPDATE genome SET pathoplexus_accession='PP_PUB1.1' WHERE accession='A1'")
    con.commit()
    con.close()
    ok = _priv_files([{"PPX_accession": "LAB_001", "strain": "s1", "date": "2024-03", "country": "Uganda"}], {"LAB_001": "ACGT" * 30})
    rep = private.check(1, *ok)
    assert rep["ok"] and rep["new"] == 1 and rep["rows"][0]["date"] == "2024-03-XX" and rep["rows"][0]["region"] == "Africa" and rep["rows"][0]["length"] == "120"
    assert any("2024-03-XX" in w for w in rep["warnings"])                       # normalisations are reported, not silent
    cases = {
        "duplicate id": (_priv_files([{"PPX_accession": "N1", "strain": "a", "date": "2024"}, {"PPX_accession": "N1", "strain": "b", "date": "2024"}], {"N1": "ACGT"}), "겹칩니다"),
        "collides with a public accession": (_priv_files([{"PPX_accession": "PP_PUB1", "strain": "a", "date": "2024"}], {"PP_PUB1": "ACGT"}), "공개 데이터"),
        "public INSDC accession": (_priv_files([{"PPX_accession": "A1", "strain": "a", "date": "2024"}], {"A1": "ACGT"}), "공개 데이터"),
        "row without a sequence": (_priv_files([{"PPX_accession": "N2", "strain": "a", "date": "2024"}], {}), "서열이 없습니다"),
        "sequence without a row": (_priv_files([{"PPX_accession": "N3", "strain": "a", "date": "2024"}], {"N3": "ACGT", "ORPHAN": "ACGT"}), "ORPHAN"),
        "not nucleotides": (_priv_files([{"PPX_accession": "N4", "strain": "a", "date": "2024"}], {"N4": "ACGTEEPP"}), "뉴클레오타이드"),
        "bad id": (_priv_files([{"PPX_accession": "a b/c", "strain": "a", "date": "2024"}], {"a": "ACGT"}), "허용되지 않는"),
    }
    for name, (files, needle) in cases.items():
        r = private.check(1, *files)
        assert not r["ok"] and any(needle in e for e in r["errors"]), (name, r["errors"])
    assert not private.check(1, "strain\tdate\nx\t2024\n", ">x\nACGT\n")["ok"]                 # no id column
    import pytest
    with pytest.raises(ValueError):
        private.import_texts(1, *cases["duplicate id"][0])                                      # errors stop the import: all or nothing
    assert private.list_rows(1) == []


def test_private_store_is_separate_and_exported_only_on_request(tmp_path, monkeypatch):
    import json
    from genomecurator import export, private, qc
    _qc_env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    con = connect()
    con.execute("UPDATE genome SET pathoplexus_accession='PP_A1.1' WHERE accession='A1'")
    con.execute("UPDATE genome SET pathoplexus_accession='PP_X2.1' WHERE accession='X2'")
    con.commit()
    con.close()
    qc.dataset_update("nextstrain/mpox/all-clades")
    files = _priv_files([{"PPX_accession": "LAB_001", "strain": "s1", "date": "2024-03-05", "country": "Uganda", "division": "Kampala", "host": "Homo sapiens"},
                         {"PPX_accession": "LAB_002", "strain": "s2", "date": "2024", "country": "Germany", "clade": "IIb"}],
                        {"LAB_001": "ACGT" * 50, "LAB_002": "ACGT" * 40 + "N" * 8})
    r0 = private.import_texts(1, *files, source_name="batch1.tsv")
    assert (r0["new"], r0["update"], r0["unchanged"]) == (2, 0, 2 - 2) and r0["versions"]["LAB_001"] == {"status": "new", "version": 1}
    r1 = private.import_texts(1, *files)
    assert (r1["new"], r1["update"], r1["unchanged"]) == (0, 0, 2)                              # the same content again: nothing changes, no new version
    # nothing leaks into the public side
    con = connect()
    assert con.execute("SELECT COUNT(*) FROM genome WHERE accession LIKE 'LAB_%'").fetchone()[0] == 0
    con.close()
    assert [x["PPX_accession"] for x in export.table(1, "nextstrain")[1]] == ["PP_A1", "PP_X2"]
    qc.run(1)
    assert all(not k.startswith("private:") for k in qc.per_genome(1)) and qc.summary(1)["by"]["seqstats"]["good"] + qc.summary(1)["by"]["seqstats"].get("bad", 0) == 2
    r = qc.run_private(1)
    assert r["stats"] == 2 and r["nextclade"] == 2
    rows = {x["id"]: x for x in private.list_rows(1)}
    assert rows["LAB_001"]["qc_nextclade"] == "good" and rows["LAB_001"]["clade"] == "Ib" and not rows["LAB_001"]["clade_differs"]
    assert rows["LAB_002"]["clade_differs"]                                                      # supplied IIb, Nextclade says Ib
    # export: public files are private-free; private files only on request, separate, wired to the workflow's keys
    out = tmp_path / "pub"
    r1 = export.write(1, out, "nextstrain")
    assert "private" not in " ".join(r1["files"]) and not (out / "metadata-private.tsv").exists() and "private_metadata" not in (out / "nextstrain-config.yaml").read_text(encoding="utf-8")
    assert json.loads((out / "manifest.json").read_text(encoding="utf-8"))["contains_private"] is False
    out2 = tmp_path / "withpriv"
    r2 = export.write(1, out2, "nextstrain", include_private=True)
    assert r2["private"] == 2 and "metadata-private.tsv" in r2["files"]
    assert "LAB_" not in (out2 / "metadata.tsv").read_text(encoding="utf-8") and "LAB_" not in (out2 / "sequences.fasta").read_text(encoding="utf-8")
    lines = (out2 / "metadata-private.tsv").read_text(encoding="utf-8").splitlines()
    cols = lines[0].split("\t")
    assert cols[:len(export.NEXTSTRAIN_COLUMNS)] == export.NEXTSTRAIN_COLUMNS
    prow = [dict(zip(cols, ln.split("\t"))) for ln in lines[1:]]
    assert [p["accession"] for p in prow] == [p["PPX_accession"] for p in prow] == ["LAB_001", "LAB_002"]    # the column the workflow's merge script matches on
    assert prow[0]["date"] == "2024-03-05" and prow[0]["region"] == "Africa" and prow[0]["dataUseTerms"] == "PRIVATE" and prow[0]["QC_rare_mutations"] == "good"
    assert prow[1]["clade"] == "IIb" and prow[1]["date"] == "2024-XX-XX"                         # the supplied clade stays; the mismatch is shown in the private tab
    assert [ln[1:] for ln in (out2 / "sequences-private.fasta").read_text(encoding="utf-8").splitlines() if ln.startswith(">")] == ["LAB_001", "LAB_002"]
    cfg = (out2 / "nextstrain-config.yaml").read_text(encoding="utf-8")
    assert 'private_metadata: "' in cfg and 'private_sequences: "' in cfg and "use_remote_files: false" in cfg
    assert "*private*" in (out2 / ".gitignore").read_text(encoding="utf-8") and (out2 / "PRIVATE_DATA.txt").exists()
    mf = json.loads((out2 / "manifest.json").read_text(encoding="utf-8"))
    assert mf["contains_private"] and mf["private"]["rows"] == 2 and mf["private"]["files"]["metadata-private.tsv"]["sha256"]
    assert export.preview(1, "nextstrain", include_private=True)["private"] == 2 and export.preview(1, "nextstrain")["private"] == 0
    assert private.delete(1, "LAB_001") and [x["id"] for x in private.list_rows(1)] == ["LAB_002"] and not private.delete(1, "LAB_001")
    import sqlite3
    c = sqlite3.connect(str(tmp_path / "g.db"))
    assert c.execute("SELECT COUNT(*) FROM qc_result WHERE accession='private:LAB_001'").fetchone()[0] == 0
    assert "private_genome" not in __import__("genomecurator.backup", fromlist=["x"]).KEEP_TABLES              # not in a general backup


def test_private_api(tmp_path, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from genomecurator import web
    _qc_env(tmp_path, monkeypatch)
    c = TestClient(web.create_app(scheduler=False))
    mt, fa = _priv_files([{"PPX_accession": "LAB_9", "strain": "s", "date": "2024-01-02", "country": "Brazil"}], {"LAB_9": "ACGT" * 20})
    bad = _priv_files([{"PPX_accession": "LAB_9", "strain": "s", "date": "2024"}], {})
    assert c.get("/private").status_code == 200
    r = c.post("/api/private/1/check", files={"metadata": ("m.tsv", bad[0]), "fasta": ("f.fa", bad[1])}).json()
    assert r["ok"] is False and r["errors"]
    r = c.post("/api/private/1/check", files={"metadata": ("m.tsv", mt), "fasta": ("f.fa", fa)}).json()
    assert r["ok"] and r["new"] == 1
    assert c.get("/api/private/1").json() == []                                                 # a check stores nothing
    assert c.post("/api/private/1/import", files={"metadata": ("m.tsv", bad[0]), "fasta": ("f.fa", bad[1])}).status_code == 400
    assert c.post("/api/private/1/import", files={"metadata": ("m.tsv", mt), "fasta": ("f.fa", fa)}).json()["new"] == 1
    assert [x["id"] for x in c.get("/api/private/1").json()] == ["LAB_9"] and "LAB_9" in c.get("/private?watch=1").text
    assert c.post("/api/private/1/check").status_code == 400 and c.post("/api/private/99/check", files={"metadata": ("m", mt), "fasta": ("f", fa)}).status_code == 404
    jid = c.post("/api/private/1/qc").json()["job"]
    for _ in range(100):
        j = c.get(f"/api/job/{jid}").json()
        if j["done"]:
            break
        time.sleep(0.05)
    assert j["status"] == "done" and j["result"]["stats"] == 1
    assert c.get("/api/export/1/preview?include_private=1").json()["private"] == 1 and c.get("/api/export/1/preview").json()["private"] == 0
    assert c.delete("/api/private/1/LAB_9").json() == {"deleted": True}


def test_backup_keeps_private_data_out_unless_asked(tmp_path, monkeypatch):
    import zipfile
    from genomecurator import backup, config, private
    src = tmp_path / "src"
    src.mkdir()
    _qc_env(src, monkeypatch)
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_B1", "strain": "s", "date": "2024-01-02", "country": "Uganda"}], {"LAB_B1": "ACGT" * 25}))
    plain = backup.export_backup(tmp_path / "plain.zip")
    assert plain["contains_private"] is False and "private_genome" not in plain["tables"]
    assert not [n for n in zipfile.ZipFile(tmp_path / "plain.zip").namelist() if n.startswith("private/")]
    full = backup.export_backup(tmp_path / "full.zip", with_private=True)
    assert full["contains_private"] is True and full["tables"]["private_genome"] == 1 and full["tables"]["private_version"] == 1 and full["private_files"] == 2     # current file + its version copy
    dst = tmp_path / "dst"
    dst.mkdir()
    monkeypatch.setattr(config, "DATA_DIR", dst)
    monkeypatch.setattr(config, "DB_PATH", dst / "g.db")
    backup.restore_backup(tmp_path / "full.zip")
    assert [r["id"] for r in private.list_rows(1)] == ["LAB_B1"] and (dst / "private" / "Mpox" / "LAB_B1.fna").exists()
    assert [h["version"] for h in private.history(1, "LAB_B1")] == [1] and (dst / "private" / "Mpox" / "_versions" / "LAB_B1.v1.fna").exists()


def test_check_drops_genomes_below_the_coverage_floor(tmp_path, monkeypatch):
    import json
    from genomecurator import config, curate, sources
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    monkeypatch.setattr(config, "MIN_COVERAGE_PCT", 80.0)
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'V','Some virus',5)")
    con.execute("INSERT INTO genome(genome_id,watch_id,accession,genome_length,genome_status) VALUES('g1',1,'OLDSHORT',600,'Partial')")     # held before, 60% of the reference
    con.execute("INSERT INTO genome(genome_id,watch_id,accession,genome_length,genome_status) VALUES('g2',1,'OLDGOOD',950,'Partial')")
    con.commit()
    con.close()
    rec = lambda ln, st="Partial", ref="": {"genome_length": ln, "genome_status": st, "refseq_accessions": ref, "accession_version": 1}
    ncbi = {"REF": rec(1000, "Complete", "NC_1"), "OKAY": rec(850), "SHORTNEW": rec(700), "TINY": rec(10), "SHORTCOMPLETE": rec(700, "Complete"),
            "OLDSHORT": rec(600), "OLDGOOD": rec(950)}
    monkeypatch.setattr(sources, "ncbi_species_scope", lambda tid: None)
    monkeypatch.setattr(sources, "pathoplexus_resolve", lambda *a, **k: None)
    monkeypatch.setattr(sources, "fetch_all", lambda *a, **k: ({"ncbi": json.loads(json.dumps(ncbi)), "bvbrc": {}, "pathoplexus": {}}, []))
    cs = curate.check("V")
    con = connect()
    items = {r["accession"]: dict(r) for r in con.execute("SELECT * FROM changeset_item WHERE changeset_id=?", (cs,))}
    summ = json.loads(con.execute("SELECT summary_json FROM changeset WHERE id=?", (cs,)).fetchone()[0])
    con.close()
    assert "SHORTNEW" not in items and "TINY" not in items                                   # never proposed
    assert {"REF", "OKAY", "SHORTCOMPLETE"} <= set(items)                                    # NCBI 'Complete' is exempt, 85% stays
    r = items["OLDSHORT"]
    assert r["kind"] == "removed" and r["decision"] == "approved" and "60.0%" in json.loads(r["provenance_json"])["reason"]   # held ones: removal pre-approved, with the reason
    assert "OLDGOOD" not in items or items["OLDGOOD"]["kind"] != "removed"
    assert summ["below_coverage"] == 3 and any("80%" in n for n in summ["notes"])
    # off switch
    monkeypatch.setattr(config, "MIN_COVERAGE_PCT", 0.0)
    assert curate._below_coverage(connect(), 1, {"ncbi": ncbi}) == {}


def _xlsx(rows, cols=("PPX_accession", "strain", "date", "country", "clade")):
    import io
    from openpyxl import Workbook
    wb = Workbook()
    ws = wb.active
    ws.title = "genomes"
    ws.append(list(cols))
    for r in rows:
        ws.append([r.get(c, "") for c in cols])
    ws.append([None] * len(cols))                       # an empty row at the end is ignored
    wb.create_sheet("예시").append(["not", "read"])
    b = io.BytesIO()
    wb.save(b)
    return b.getvalue()


def test_private_spreadsheet_template_and_xlsx_upload(tmp_path, monkeypatch):
    import datetime
    from genomecurator import private
    _qc_env(tmp_path, monkeypatch)
    from openpyxl import load_workbook
    import io
    wb = load_workbook(io.BytesIO(private.template_xlsx("Mpox", "Orthopoxvirus monkeypox")))
    assert wb.sheetnames == ["genomes", "열 설명", "예시"]
    assert [c.value for c in wb["genomes"][1]] == [c[0] for c in private.COLUMNS]
    assert wb["genomes"].data_validations.dataValidation and "Ib" in wb["genomes"].data_validations.dataValidation[0].formula1      # clade drop-down for mpox
    assert [r[0].value for r in wb["열 설명"].iter_rows(min_row=2, max_row=4)] == ["PPX_accession", "strain", "date"]
    # a filled spreadsheet: real Excel dates become ISO text, numbers stay readable, the example / explanation sheets are not read
    data = _xlsx([{"PPX_accession": "XL_1", "strain": "s1", "date": datetime.date(2025, 3, 14), "country": "Uganda"},
                  {"PPX_accession": "XL_2", "strain": "s2", "date": "2025-04", "country": "Kenya", "clade": "IIb"}])
    fa = ">XL_1\n" + "ACGT" * 30 + "\n>XL_2\n" + "ACGT" * 20 + "\n"
    rep = private.check(1, data, fa, defaults={"clade": "Ib"}, filename="batch.xlsx")
    assert rep["ok"] and rep["new"] == 2
    rows = {r["_id"]: r for r in rep["rows"]}
    assert rows["XL_1"]["date"] == "2025-03-14" and rows["XL_1"]["clade"] == "Ib"             # blank clade -> the batch default
    assert rows["XL_2"]["date"] == "2025-04-XX" and rows["XL_2"]["clade"] == "IIb"            # a supplied clade is never overwritten by the default
    assert private.import_texts(1, data, fa, "batch.xlsx", defaults={"clade": "Ib"})["new"] == 2
    csvtext = "PPX_accession,strain,date,country\nCS_1,s,2025-01-02,Peru\n"                  # comma separated text works too
    assert private.check(1, csvtext, ">CS_1\n" + "ACGT" * 10 + "\n")["ok"]


def test_private_one_by_one_and_versions(tmp_path, monkeypatch):
    import pytest
    from genomecurator import private
    _qc_env(tmp_path, monkeypatch)
    fields = {"PPX_accession": "ONE_1", "strain": "a", "date": "2025-05-06", "country": "Uganda", "division": "Kampala"}
    seq = "ACGT" * 40
    assert private.one(1, fields, seq, dry_run=True)["new"] == 1 and private.list_rows(1) == []     # a check stores nothing
    r = private.one(1, fields, ">ONE_1 pasted\n" + seq + "\n", note="first")                         # FASTA text or a bare sequence
    assert r["new"] == 1 and r["versions"]["ONE_1"] == {"status": "new", "version": 1}
    assert private.one(1, fields, seq)["unchanged"] == 1 and len(private.history(1, "ONE_1")) == 1  # same content: no new version
    r2 = private.one(1, {**fields, "division": "Wakiso", "host": "Homo sapiens"}, seq, note="fix place")
    assert r2["update"] == 1 and r2["versions"]["ONE_1"]["version"] == 2
    r3 = private.one(1, {**fields, "division": "Wakiso", "host": "Homo sapiens"}, "ACGT" * 41, note="resequenced")
    assert r3["versions"]["ONE_1"]["version"] == 3
    h = private.history(1, "ONE_1")
    assert [x["version"] for x in h] == [3, 2, 1] and [x["current"] for x in h] == [True, False, False]
    assert h[1]["changes"] == ["division", "host"] and not h[1]["sequence_changed"] and h[1]["note"] == "fix place"
    assert h[0]["sequence_changed"] and h[2]["new"]
    row = private.list_rows(1)[0]
    assert row["version"] == 3 and row["length"] == 164
    files = sorted(p.name for p in (tmp_path / "private" / "Mpox" / "_versions").iterdir())
    assert files == ["ONE_1.v1.fna", "ONE_1.v3.fna"]                                                 # v2 shares v1's sequence, v3 has its own
    # restore: the old version comes back as a NEW version, history is never rewritten
    rr = private.restore(1, "ONE_1", 1)
    assert rr["version"] == 4 and private.list_rows(1)[0]["length"] == 160
    assert [x["version"] for x in private.history(1, "ONE_1")][:2] == [4, 3] and private.history(1, "ONE_1")[0]["note"] == "버전 1 복원"
    # delete keeps the history; the genome can be brought back
    assert private.delete(1, "ONE_1") and private.list_rows(1) == []
    hh = private.history(1, "ONE_1")
    assert hh[0]["deleted"] and hh[0]["version"] == 5
    assert private.restore(1, "ONE_1", 3)["version"] == 6 and private.list_rows(1)[0]["length"] == 164
    with pytest.raises(ValueError):
        private.one(1, fields, ">a\nACGT\n>b\nACGT\n")                                              # one sequence only
    with pytest.raises(ValueError):
        private.restore(1, "ONE_1", 99)


def test_private_import_respects_the_coverage_floor(tmp_path, monkeypatch):
    from genomecurator import private
    _qc_env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    con = connect()
    con.execute("UPDATE genome SET refseq_accessions='NC_1', genome_length=1000 WHERE accession='A1'")      # a reference of 1,000 bp
    con.commit()
    con.close()
    short = private.one(1, {"PPX_accession": "SHORT_1", "strain": "s", "date": "2025", "country": "Uganda"}, "ACGT" * 150, dry_run=True)   # 600 bp = 60%
    assert not short["ok"] and any("커버리지 하한" in e and "60.0%" in e for e in short["errors"])
    good = private.one(1, {"PPX_accession": "GOOD_1", "strain": "s", "date": "2025", "country": "Uganda"}, "ACGT" * 230, dry_run=True)      # 920 bp
    assert good["ok"]
    monkeypatch.setattr(__import__("genomecurator.config", fromlist=["x"]), "MIN_COVERAGE_PCT", 0.0)
    assert private.one(1, {"PPX_accession": "SHORT_1", "strain": "s", "date": "2025", "country": "Uganda"}, "ACGT" * 150, dry_run=True)["ok"]


def test_private_parse_endpoint_opens_a_file_in_the_sheet(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from genomecurator import web
    _qc_env(tmp_path, monkeypatch)
    c = TestClient(web.create_app(scheduler=False))
    x = _xlsx([{"PPX_accession": "P1", "strain": "s", "date": "2025-01-02", "country": "Peru"}])
    r = c.post("/api/private/parse", files={"metadata": ("a.xlsx", x)}).json()
    assert r["cols"] == ["PPX_accession", "strain", "date", "country", "clade"] and r["rows"][0]["strain"] == "s" and len(r["rows"]) == 1
    r = c.post("/api/private/parse", files={"metadata": ("a.tsv", b"PPX_accession\tstrain\nQ1\tz\n")}).json()
    assert r["rows"] == [{"PPX_accession": "Q1", "strain": "z"}]
    assert c.get("/private?watch=1").text.count("jspreadsheet") >= 2 and c.post("/api/private/parse").status_code == 400


def test_private_api_one_template_and_history(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from genomecurator import web
    _qc_env(tmp_path, monkeypatch)
    c = TestClient(web.create_app(scheduler=False))
    t = c.get("/api/private/template.xlsx?watch=1")
    assert t.status_code == 200 and t.content[:2] == b"PK" and "spreadsheetml" in t.headers["content-type"]
    assert c.get("/api/private/template.tsv").text.startswith("PPX_accession\tstrain\tdate")
    form = {"PPX_accession": "API_1", "strain": "s", "date": "2025-02", "country": "Brazil", "sequence": "ACGT" * 30, "note": "via form", "default_clade": "Ib"}
    r = c.post("/api/private/1/one?dry_run=1", data=form).json()
    assert r["ok"] and r["new"] == 1 and c.get("/api/private/1").json() == []
    r = c.post("/api/private/1/one", data=form).json()
    assert r["new"] == 1 and r["versions"]["API_1"]["version"] == 1
    assert c.get("/api/private/1").json()[0]["clade"] == "Ib" and c.get("/api/private/1").json()[0]["version"] == 1
    r = c.post("/api/private/1/one", data={**form, "country": "Peru"}).json()
    assert r["update"] == 1
    h = c.get("/api/private/1/history?id=API_1").json()
    assert [x["version"] for x in h] == [2, 1] and h[0]["changes"] == ["country"]
    assert c.post("/api/private/1/API_1/restore?version=1").json()["version"] == 3
    assert c.post("/api/private/1/API_1/restore?version=99").status_code == 400
    assert c.post("/api/private/1/one", data={"strain": "x", "sequence": "ACGT"}).status_code == 400                # no id
    # a spreadsheet through the API
    data = _xlsx([{"PPX_accession": "API_2", "strain": "s", "date": "2025", "country": "Peru"}])
    r = c.post("/api/private/1/check", files={"metadata": ("b.xlsx", data), "fasta": ("f.fa", ">API_2\n" + "ACGT" * 20 + "\n")}, data={"default_clade": "Ib"}).json()
    assert r["ok"] and r["new"] == 1
    assert "비공개" in c.get("/private?watch=1").text and "API_1" in c.get("/private?watch=1").text


def test_exclusion_is_a_mark_not_a_deletion(tmp_path, monkeypatch):
    import json
    import pytest
    from genomecurator import exclusions, export, private, qc
    _qc_env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    con = connect()
    con.execute("UPDATE genome SET pathoplexus_accession='PP_A1.1' WHERE accession='A1'")
    con.execute("UPDATE genome SET pathoplexus_accession='PP_X2.1' WHERE accession='X2'")
    con.commit()
    con.close()
    qc.dataset_update("nextstrain/mpox/all-clades")
    qc.run(1)
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_E1", "strain": "a", "date": "2025", "country": "Uganda"},
                                          {"PPX_accession": "LAB_E2", "strain": "b", "date": "2025", "country": "Uganda"}], {"LAB_E1": "ACGT" * 40, "LAB_E2": "ACGT" * 41}))
    assert [x["PPX_accession"] for x in export.table(1, "nextstrain")[1]] == ["PP_A1", "PP_X2"]
    r = exclusions.add(1, ["X2", "NOPE", "X2"], "public", "bad assembly")
    assert r == {"excluded": 1, "new": 1, "unknown": ["NOPE"]}                                # unknown ids are reported, never stored
    assert exclusions.add(1, ["X2"], "public", "still bad")["new"] == 0                       # again: only the reason is updated
    assert exclusions.marked(1)["X2"]["reason"] == "still bad"
    with pytest.raises(ValueError):
        exclusions.add(1, ["X2"], "elsewhere")
    # nothing was deleted: the genome, its sequence, its QC are all still there
    con = connect()
    assert con.execute("SELECT COUNT(*) FROM genome WHERE accession='X2'").fetchone()[0] == 1
    assert con.execute("SELECT file_path FROM genome_sequence WHERE accession='X2'").fetchone()[0]
    assert con.execute("SELECT COUNT(*) FROM qc_result WHERE accession='X2'").fetchone()[0] == 2
    con.close()
    assert (tmp_path / "genomes" / "X2.fna").exists()
    # ...but the dataset leaves it out, says so, and the preview warns
    assert [x["PPX_accession"] for x in export.table(1, "nextstrain")[1]] == ["PP_A1"]
    assert export.preview(1, "nextstrain")["excluded"] == 1 and any("제외하도록 표시" in w for w in export.preview(1, "nextstrain")["warnings"])
    out = tmp_path / "ds"
    res = export.write(1, out, "nextstrain")
    assert res["rows"] == 1 and "excluded.tsv" in res["files"]
    lines = (out / "excluded.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[0] == "accession\tstrain\treason\texcluded_at" and lines[1].split("\t")[:3] == ["X2", "", "still bad"]
    mf = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert mf["excluded"]["public"] == 1 and mf["excluded"]["files"] == ["excluded.tsv"]
    assert "PP_X2" not in (out / "metadata.tsv").read_text(encoding="utf-8") and "PP_X2" not in (out / "sequences.fasta").read_text(encoding="utf-8")
    # the private store too (and its file only travels with the private export)
    assert exclusions.add(1, ["LAB_E2"], "private", "contaminated")["excluded"] == 1 and exclusions.add(1, ["PP_A1"], "private")["unknown"] == ["PP_A1"]
    assert [x["id"] for x in private.list_rows(1) if x["excluded"]] == ["LAB_E2"]
    assert [r["PPX_accession"] for r in export.private_table(1)[1]] == ["LAB_E1"]
    out2 = tmp_path / "ds2"
    r2 = export.write(1, out2, "nextstrain", include_private=True)
    assert r2["private"] == 1 and "excluded-private.tsv" in r2["files"] and "LAB_E2" not in (out2 / "metadata-private.tsv").read_text(encoding="utf-8")
    assert not (out / "excluded-private.tsv").exists()
    assert len(private.list_rows(1)) == 2 and (tmp_path / "private" / "Mpox" / "LAB_E2.fna").exists()      # still stored
    # lifting the mark brings them back, nothing else changes
    assert exclusions.remove(1, ["X2"]) == {"included": 1} and exclusions.remove(1, ["X2"]) == {"included": 0}
    assert [x["PPX_accession"] for x in export.table(1, "nextstrain")[1]] == ["PP_A1", "PP_X2"]
    assert export.table(1, "nextstrain", skip_excluded=True)[1] and export.write(1, tmp_path / "ds3", "nextstrain")["excluded"] == 0
    assert not (tmp_path / "ds3" / "excluded.tsv").exists()


def test_exclusion_api_and_pages(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from genomecurator import private, web
    _qc_env(tmp_path, monkeypatch)
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_Z1", "strain": "a", "date": "2025", "country": "Uganda"}], {"LAB_Z1": "ACGT" * 40}))
    c = TestClient(web.create_app(scheduler=False))
    r = c.post("/api/exclusions/1", json={"ids": ["A1", "ZZ"], "scope": "public", "action": "exclude", "reason": "short read"}).json()
    assert r["excluded"] == 1 and r["unknown"] == ["ZZ"]
    assert [x["id"] for x in c.get("/api/exclusions/1").json()] == ["A1"] and c.get("/api/exclusions/1?scope=private").json() == []
    page = c.get("/genomes?watch=1").text
    import json, re
    data = json.loads(re.search(r'<script type="application/json" id="gdata">(.*?)</script>', page, re.S).group(1))
    assert "exclSel" in page and {r["accession"]: r["excluded"] for r in data["rows"]} == {"A1": "제외", "X2": "포함"}
    assert {r["accession"]: r["excl_reason"] for r in data["rows"]}["A1"] == "short read" and any(f["k"] == "excluded" for f in data["facets"])
    assert c.post("/api/exclusions/1", json={"ids": ["LAB_Z1"], "scope": "private", "action": "exclude"}).json()["excluded"] == 1
    assert "제외됨" in c.get("/private?watch=1").text
    assert c.post("/api/exclusions/1", json={"ids": ["A1"], "scope": "public", "action": "include"}).json() == {"included": 1}
    assert c.post("/api/exclusions/99", json={"ids": ["A1"], "scope": "public", "action": "exclude"}).status_code == 400
    assert c.get("/api/export/1/preview?profile=nextstrain").json()["excluded"] == 0


def _ds_env(tmp_path, monkeypatch):
    """_qc_env plus: A1 = clade Ib (2024, Uganda), X2 = clade IIa (2019, Peru), N3 = clade Ib without a sequence file; all with PPX accessions; QC run."""
    from genomecurator import qc
    _qc_env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    con = connect()
    con.execute("UPDATE genome SET pathoplexus_accession='PP_A1.1', clade='Ib', isolation_country='Uganda', collection_date='2024-05-06', accession_version=1 WHERE accession='A1'")
    con.execute("UPDATE genome SET pathoplexus_accession='PP_X2.1', clade='IIa', isolation_country='Peru', collection_date='2019-01-02', accession_version=1 WHERE accession='X2'")
    con.execute("INSERT INTO genome(genome_id,watch_id,accession,genome_length,clade,pathoplexus_accession,isolation_country,collection_date,accession_version) VALUES('g9',1,'N3',200,'Ib','PP_N3.1','Uganda','2023',1)")
    con.commit()
    con.close()
    qc.dataset_update("nextstrain/mpox/all-clades")
    qc.run(1)


def test_dataset_rules_and_storage(tmp_path, monkeypatch):
    import pytest
    from genomecurator import datasets, export
    _ds_env(tmp_path, monkeypatch)
    assert datasets.describe({}) == "제한 없음 (전체)" and "clade Ib" in datasets.describe({"clade": ["Ib"]}) and "서열 없는 것 포함" in datasets.describe({"require_sequence": False})
    n = datasets.normalize({"clade": [" Ib ", ""], "year_min": "2020", "year_max": "", "bogus": 1, "min_length": "abc"})
    assert n["clade"] == ["Ib"] and n["year_min"] == 2020 and n["year_max"] is None and n["min_length"] is None and "bogus" not in n and n["require_sequence"] is True
    ib = datasets.save(1, "clade Ib", {"clade": ["Ib"]}, "only Ib", build="clade-i")
    assert ib["id"] and ib["rules"]["clade"] == ["Ib"] and ib["build"] == "clade-i" and datasets.find(1, "clade Ib")["id"] == ib["id"]
    with pytest.raises(ValueError):
        datasets.save(1, "clade Ib", {})                                                         # names are unique per virus
    with pytest.raises(ValueError):
        datasets.save(1, "x", {}, build="nope")
    with pytest.raises(ValueError):
        datasets.save(1, "  ", {})
    ids = lambda ds, **kw: [r["PPX_accession"] for r in export.table(1, "nextstrain", ds, **kw)[1]]
    assert ids(ib) == ["PP_A1"]                                                                  # N3 is Ib too, but has no sequence file
    both = datasets.save(1, "Ib incl. no sequence", {"clade": ["Ib"], "require_sequence": False})
    assert ids(both) == ["PP_A1", "PP_N3"]
    assert ids(datasets.save(1, "years", {"year_min": 2020})) == ["PP_A1"] and ids(datasets.save(1, "old", {"year_max": 2019})) == ["PP_X2"]
    assert ids(datasets.save(1, "peru", {"country": ["Peru"]})) == ["PP_X2"] and ids(datasets.save(1, "len", {"min_length": 190})) == ["PP_A1", "PP_X2"] and ids(datasets.save(1, "len2", {"min_length": 201})) == []
    assert ids(datasets.save(1, "all", {})) == ["PP_A1", "PP_X2"]
    # manual add / remove: add wins over the rules, remove wins over everything; ids may be the PPX or the INSDC accession
    assert ids(datasets.save(1, "Ib + X2", {"clade": ["Ib"]}, manual_include=["X2"])) == ["PP_A1", "PP_X2"]
    assert ids(datasets.save(1, "Ib - A1", {"clade": ["Ib"]}, manual_exclude=["PP_A1"])) == []
    # a QC status can leave members out (X2's Nextclade QC is bad in the fake)
    assert ids(datasets.save(1, "all but bad", {"qc_not": ["bad"]})) == ["PP_A1"]
    # the global exclusion still applies to every dataset
    from genomecurator import exclusions
    exclusions.add(1, ["A1"], "public", "x")
    assert ids(ib) == [] and ids(datasets.save(1, "all2", {})) == ["PP_X2"]
    exclusions.remove(1, ["A1"])
    # update in place; the definition can be deleted, the genomes stay
    up = datasets.save(1, "clade Ib (renamed)", {"clade": ["Ib", "IIa"]}, dataset_id=ib["id"])
    assert up["id"] == ib["id"] and datasets.get(ib["id"])["name"] == "clade Ib (renamed)" and ids(up) == ["PP_A1", "PP_X2"]
    assert datasets.delete(ib["id"]) and not datasets.delete(ib["id"]) and datasets.get(ib["id"]) is None
    from genomecurator.db import connect
    assert connect().execute("SELECT COUNT(*) FROM genome WHERE accession='A1'").fetchone()[0] == 1


def test_dataset_export_snapshot_and_private(tmp_path, monkeypatch):
    import json
    from genomecurator import datasets, export, private
    _ds_env(tmp_path, monkeypatch)
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_U1", "strain": "a", "date": "2025", "country": "Uganda"},
                                          {"PPX_accession": "LAB_P1", "strain": "b", "date": "2025", "country": "Peru"}], {"LAB_U1": "ACGT" * 40, "LAB_P1": "ACGT" * 41}))
    ds = datasets.save(1, "Uganda only", {"country": ["Uganda"]}, "test", include_private=True, build="clade-i")
    pv = export.preview(1, "nextstrain", include_private=True, dataset=ds)
    assert pv["selected"] == 1 and pv["private"] == 1                                            # the rules apply to private rows too
    r = export.write(1, None, "nextstrain", dataset=ds)
    out = __import__("pathlib").Path(r["dir"])
    assert out.parent.name == "Uganda_only" and out.parent.parent.name == "Mpox"                 # exports/<virus>/<dataset>/<stamp>
    assert r["rows"] == 1 and r["private"] == 1                                                  # the dataset asked for the private part, no flag needed
    assert "PP_X2" not in (out / "metadata.tsv").read_text(encoding="utf-8") and "LAB_U1" in (out / "metadata-private.tsv").read_text(encoding="utf-8")
    assert "LAB_P1" not in (out / "metadata-private.tsv").read_text(encoding="utf-8")
    mf = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert mf["dataset"]["name"] == "Uganda only" and mf["dataset"]["rules"]["country"] == ["Uganda"] and mf["dataset"]["members_sha256"] == datasets.members_hash(["PP_A1", "private:LAB_U1"])
    assert "defaults/clade-i/config.yaml" in (out / "nextstrain-config.yaml").read_text(encoding="utf-8")
    h = datasets.export_history(ds["id"])
    assert len(h) == 1 and h[0]["rows"] == 1 and h[0]["private_rows"] == 1 and h[0]["dir"] == str(out) and h[0]["members_sha256"] == mf["dataset"]["members_sha256"]
    assert datasets.listing(1)[0]["last_export"] and export.list_exports()[0]["dataset"] == "Uganda only"
    # the same dataset exported again later has the same member hash while the members are the same (traceable), and differs when they change
    r2 = export.write(1, tmp_path / "again", "nextstrain", dataset=ds)
    assert json.loads((tmp_path / "again" / "manifest.json").read_text(encoding="utf-8"))["dataset"]["members_sha256"] == mf["dataset"]["members_sha256"] and r2["rows"] == 1
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_U2", "strain": "c", "date": "2025", "country": "Uganda"}], {"LAB_U2": "ACGT" * 42}))
    export.write(1, tmp_path / "later", "nextstrain", dataset=ds)
    assert json.loads((tmp_path / "later" / "manifest.json").read_text(encoding="utf-8"))["dataset"]["members_sha256"] != mf["dataset"]["members_sha256"]
    # a plain (non-dataset) export is untouched
    assert export.write(1, tmp_path / "plain", "nextstrain")["rows"] == 2 and json.loads((tmp_path / "plain" / "manifest.json").read_text(encoding="utf-8"))["dataset"] is None


def test_dataset_api(tmp_path, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from genomecurator import web
    _ds_env(tmp_path, monkeypatch)
    c = TestClient(web.create_app(scheduler=False))
    assert c.get("/datasets?watch=1").status_code == 200
    f = c.get("/api/datasets/facets?watch=1").json()
    cl = {x["value"]: x for x in f["fields"]["clade"]}
    assert f["genomes"] == 3 and f["with_sequence"] == 2 and cl["Ib"]["n"] == 2 and cl["Ib"]["with_seq"] == 1 and cl["IIa"]["n"] == 1
    assert c.get("/api/datasets/facets?watch=99").status_code == 404
    body = {"watch": 1, "name": "Ib", "description": "d", "build": "clade-i", "rules": {"clade": ["Ib"]}, "include_private": False, "manual_include": [], "manual_exclude": []}
    pv = c.post("/api/datasets/preview", json=body).json()
    assert pv["selected"] == 1 and pv["by_clade"] == [["Ib", 1]] and pv["by_country"][0][0] == "Uganda" and "clade Ib" in pv["summary"]
    assert c.get("/api/datasets?watch=1").json() == []                                           # a preview saves nothing
    ds = c.post("/api/datasets", json=body).json()
    assert ds["id"] and c.post("/api/datasets", json=body).status_code == 400                    # duplicate name
    assert [d["name"] for d in c.get("/api/datasets?watch=1").json()] == ["Ib"]
    t = c.get(f"/api/export/1/metadata.tsv?dataset={ds['id']}")
    fa = c.get(f"/api/export/1/sequences.fasta?dataset={ds['id']}")
    assert [ln.split("\t")[0] for ln in t.text.splitlines()[1:]] == ["PP_A1"] == [ln[1:] for ln in fa.text.splitlines() if ln.startswith(">")]
    assert c.get(f"/api/export/1/metadata.tsv?dataset={ds['id']}").status_code == 200 and c.get("/api/export/1/metadata.tsv?dataset=999").status_code == 400
    jid = c.post(f"/api/datasets/{ds['id']}/export").json()["job"]
    for _ in range(100):
        j = c.get(f"/api/job/{jid}").json()
        if j["done"]:
            break
        time.sleep(0.05)
    assert j["status"] == "done" and j["result"]["rows"] == 1 and j["kind"] == "dataset"
    assert len(c.get(f"/api/datasets/{ds['id']}/history").json()) == 1
    assert c.post("/api/datasets/999/export").status_code == 404
    assert c.delete(f"/api/datasets/{ds['id']}").json() == {"deleted": True} and c.get("/api/datasets?watch=1").json() == []


def _embedded(html, ident):
    import json
    import re
    return json.loads(re.search(rf'<script type="application/json" id="{ident}">(.*?)</script>', html, re.S).group(1))


def test_dataset_is_browsable_like_the_genome_table_and_on_the_dashboard(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from genomecurator import datasets, private, sources, web
    _ds_env(tmp_path, monkeypatch)
    monkeypatch.setattr(sources, "pathoplexus_slug_for", lambda acc: None)             # the genome page looks the Pathoplexus organism up once; no network in tests
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_V1", "strain": "kr-1", "date": "2025-02-03", "country": "Uganda", "clade": "Ib"},
                                          {"PPX_accession": "LAB_V2", "strain": "pe-1", "date": "2025", "country": "Peru"}], {"LAB_V1": "ACGT" * 40, "LAB_V2": "ACGT" * 41}))
    ib = datasets.save(1, "Ib", {"clade": ["Ib"]}, include_private=True)
    pub, priv = datasets.member_ids(1, ib)
    assert pub == {"A1", "N3"} and priv == {"LAB_V1"}                       # N3 has no sequence file: still a member for browsing (an export would skip it)
    assert datasets.member_ids(1, datasets.save(1, "pub only", {"clade": ["Ib"]}))[1] == set()          # private rows only when the dataset asks
    c = TestClient(web.create_app(scheduler=False))
    page = c.get(f"/dataset/{ib['id']}")
    assert page.status_code == 200 and "Dataset View" in page.text and "이 데이터셋의 대시보드" in page.text
    d = _embedded(page.text, "gdata")
    rows = {r["accession"]: r for r in d["rows"]}
    assert set(rows) == {"A1", "N3", "LAB_V1"} and rows["LAB_V1"]["is_private"] and rows["LAB_V1"]["source_kind"] == "비공개" and rows["A1"]["source_kind"] == "공개"
    assert rows["N3"]["has_fasta"] == "No" and rows["LAB_V1"]["has_fasta"] == "Yes" and rows["LAB_V1"]["completeness"] == "비공개"
    assert [i for i, r in enumerate(d["rows"])] == [r["id"] for r in d["rows"]]                          # row ids are positions (the table's selection uses them)
    assert any(f["k"] == "source_kind" for f in d["facets"]) and any(col["k"] == "source_kind" for col in d["cols"])
    assert "서열이 없는 1건" in page.text.replace("\n", " ")
    # the ordinary genome table is unchanged: every curated genome, no private rows, no dataset header
    g = c.get("/genomes?watch=1")
    gd = _embedded(g.text, "gdata")
    assert {r["accession"] for r in gd["rows"]} == {"A1", "X2", "N3"} and "Dataset View" not in g.text and not any(col["k"] == "source_kind" for col in gd["cols"])
    assert c.get("/dataset/999", follow_redirects=False).status_code == 303
    # the dashboard over the dataset
    dash = c.get(f"/dashboard?dataset={ib['id']}")
    dd = _embedded(dash.text, "dd")
    assert dash.status_code == 200 and dd["total"] == 3 and dd["clade"]["labels"] == ["Ib"] and "데이터셋: Ib" in dash.text and f"/dataset/{ib['id']}?" in dash.text
    full = _embedded(c.get("/dashboard?watch=1").text, "dd")
    assert full["total"] == 3 and set(full["clade"]["labels"]) == {"Ib", "IIa"} and "/genomes?watch=1&" in c.get("/dashboard?watch=1").text
    only = _embedded(c.get(f"/dashboard?dataset={datasets.find(1, 'pub only')['id']}").text, "dd")
    assert only["total"] == 2                                                                                # public members only when the dataset leaves the private store out
    assert c.get("/dashboard?dataset=999&watch=1").status_code == 200                                         # an unknown dataset falls back to everything


def test_dataset_qc_runs_only_for_its_members(tmp_path, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from genomecurator import datasets, private, qc, sources, web
    _qc_env(tmp_path, monkeypatch)                                              # A1, X2 with sequences; no QC yet
    monkeypatch.setattr(sources, "pathoplexus_slug_for", lambda acc: None)
    from genomecurator.db import connect
    con = connect()
    con.execute("UPDATE genome SET clade='Ib' WHERE accession='A1'")
    con.execute("UPDATE genome SET clade='IIa' WHERE accession='X2'")
    con.commit()
    con.close()
    qc.dataset_update("nextstrain/mpox/all-clades")
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_Q1", "strain": "a", "date": "2025", "country": "Uganda", "clade": "Ib"},
                                          {"PPX_accession": "LAB_Q2", "strain": "b", "date": "2025", "country": "Uganda", "clade": "IIb"}], {"LAB_Q1": "ACGT" * 40, "LAB_Q2": "ACGT" * 41}))
    ds = datasets.save(1, "Ib", {"clade": ["Ib"]}, include_private=True)
    s0 = qc.summary_for(1, *datasets.member_ids(1, ds))
    assert s0["public"]["members"] == 1 and s0["public"]["not_run"] == 1 and s0["private"]["members"] == 1 and s0["private"]["not_run"] == 1
    res = datasets.run_qc(ds)
    assert res["public"]["sequences"] == 1 and res["public"]["stats"] == 1 and res["private"]["sequences"] == 1
    con = connect()
    done = {r[0] for r in con.execute("SELECT DISTINCT accession FROM qc_result")}
    con.close()
    assert done == {"A1", "private:LAB_Q1"}                                      # X2 and LAB_Q2 are not in the dataset: untouched
    s1 = qc.summary_for(1, *datasets.member_ids(1, ds))
    assert s1["public"]["not_run"] == 0 and s1["public"]["seqstats"] == {"good": 1} and s1["public"]["nextclade"] == {"good": 1} and s1["private"]["seqstats"] == {"good": 1}
    assert datasets.run_qc(ds)["public"]["stats"] == 0                          # nothing changed: nothing recomputed
    assert datasets.run_qc(ds, force=True)["public"]["stats"] == 1
    # without private members in the dataset the private store is left alone
    pub_only = datasets.save(1, "Ib public", {"clade": ["Ib"]})
    assert datasets.run_qc(pub_only)["private"] is None
    c = TestClient(web.create_app(scheduler=False))
    assert c.get(f"/api/datasets/{ds['id']}/qc").json()["public"]["members"] == 1 and c.get("/api/datasets/999/qc").status_code == 404
    page = c.get(f"/dataset/{ds['id']}").text
    assert "이 데이터셋 QC 실행" in page and "모두 QC 완료" in page
    jid = c.post(f"/api/datasets/{ds['id']}/qc?force=1").json()["job"]
    for _ in range(200):
        j = c.get(f"/api/job/{jid}").json()
        if j["done"]:
            break
        time.sleep(0.05)
    assert j["status"] == "done" and j["kind"] == "dsqc" and j["result"]["public"]["sequences"] == 1
    assert c.post("/api/datasets/999/qc").status_code == 404
