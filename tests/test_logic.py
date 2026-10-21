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
    from genomecurator import datasets
    datasets.save(1, "mine", {"clade": ["Ib"]}, include_private=True, manual_include=["LAB_B1", "PUB1"], manual_exclude=["PUB2"])
    from genomecurator.db import connect as _c
    _con = _c()
    _con.execute("INSERT OR REPLACE INTO dataset_exclusion(scope,watch_id,id,reason,excluded_at) VALUES('private',1,'LAB_B1','x','now'),('public',1,'PUB9','bad','now')")
    _con.commit()
    _con.close()
    plain = backup.export_backup(tmp_path / "plain.zip")
    assert plain["contains_private"] is False and "private_genome" not in plain["tables"] and "submission_package" not in plain["tables"]
    assert plain["tables"]["dataset"] == 1 and plain["tables"]["dataset_exclusion"] == 1                                    # the user's dataset decisions are in a general backup ...
    import sqlite3 as _s, json as _j
    with zipfile.ZipFile(tmp_path / "plain.zip") as _z:
        (tmp_path / "x.db").write_bytes(_z.read("genomecurator.db"))
    _d = _s.connect(str(tmp_path / "x.db"))
    assert _j.loads(_d.execute("SELECT manual_include_json FROM dataset").fetchone()[0]) == ["PUB1"]                         # ... but without the private id
    assert _d.execute("SELECT id FROM dataset_exclusion").fetchall() == [("PUB9",)] and _d.execute("SELECT include_private FROM dataset").fetchone()[0] == 0
    _d.close()
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
    assert "exclSel" in page and {r["accession"]: r["excluded"] for r in data["rows"] if not r.get("is_private")} == {"A1": "제외", "X2": "포함"}
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
    # the ordinary genome table: every curated genome and, marked, the private ones of the virus (no dataset header)
    g = c.get("/genomes?watch=1")
    gd = _embedded(g.text, "gdata")
    assert {r["accession"] for r in gd["rows"]} == {"A1", "X2", "N3", "LAB_V1", "LAB_V2"} and "Dataset View" not in g.text and any(col["k"] == "source_kind" for col in gd["cols"])
    assert {r["accession"]: r["source_kind"] for r in gd["rows"]}["LAB_V1"] == "비공개" and {r["accession"]: r["source_kind"] for r in gd["rows"]}["X2"] == "공개"
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


def _ns_env(tmp_path, monkeypatch):
    """_ds_env plus a fake Nextstrain world: a workflow checkout (Snakefile + a clade-i build config) and a fake environment whose `nextstrain build`
    writes <auspice_name>.json the way the real workflow would, printing snakemake-style progress."""
    import os
    import stat
    from genomecurator import nextstrain, sources
    _ds_env(tmp_path, monkeypatch)
    monkeypatch.setattr(sources, "pathoplexus_slug_for", lambda acc: None)
    co = tmp_path / "wf"
    (co / "phylogenetic" / "defaults" / "clade-i").mkdir(parents=True)
    (co / "phylogenetic" / "Snakefile").write_text("")
    (co / "phylogenetic" / "defaults" / "clade-i" / "config.yaml").write_text("filter:\n  min_length: 170000\n  min_date: 1900\ntimetree: false\nroot: best\n")
    monkeypatch.setenv("GC_NEXTSTRAIN_DIR", str(co))
    monkeypatch.setenv("GC_AUSPICE_PORT", "49173")                       # never the real viewer's port: a viewer that is running on this machine must not change what the tests see
    env = tmp_path / "envs" / "gc-nextstrain"
    (env / "bin").mkdir(parents=True)
    fake = env / "bin" / "nextstrain"
    fake.write_text("""#!/usr/bin/env python3
import os, re, sys, json
a = sys.argv[1:]
if a[0] == "--version":
    print("Nextstrain CLI 9.9"); sys.exit(0)
cfgs = a[a.index("--configfile") + 1:]
cfgs = cfgs[:cfgs.index("--cores")] if "--cores" in cfgs else cfgs
overlay = open(cfgs[-1]).read()
name = re.search(r'auspice_name: "([^"]+)"', overlay).group(1)
if os.environ.get("FAKE_NS_FAIL"):
    print("rule align:"); print("Error in rule align: boom"); sys.exit(3)
for i, r in enumerate(["align", "tree", "refine", "export"], 1):
    print(f"rule {r}:"); print(f"{i} of 4 steps ({25 * i}%) done")
os.makedirs("auspice", exist_ok=True)
json.dump({"meta": {"title": name}, "tree": {"name": "root"}}, open(f"auspice/{name}.json", "w"))
json.dump({"nuc": "ACGT"}, open(f"auspice/{name}_root-sequence.json", "w"))
""")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    for t in ("auspice", "augur", "snakemake", "nextclade", "iqtree"):
        f = env / "bin" / t
        f.write_text(f"#!/bin/sh\necho '{t} 1.0'\n")
        f.chmod(f.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setattr(nextstrain, "env_prefix", lambda: env)
    return nextstrain


import os as _os
import pytest as _pytest
_unix_only = _pytest.mark.skipif(_os.name == "nt", reason="the fake `nextstrain` is a shell-style script")


def test_nextstrain_settings_validation_and_merge(tmp_path, monkeypatch):
    import pytest
    from genomecurator import datasets
    ns = _ns_env(tmp_path, monkeypatch)
    assert ns.parse_params("") == {} and ns.parse_params("# only a comment\n") == {}
    assert ns.parse_params("filter:\n  min_length: 150000\n")["filter"]["min_length"] == 150000
    for bad, why in (("inputs: []\n", "관리"), ("build_name: x\n", "관리"), ("- a\n- b\n", "YAML"), ("filter: [unclosed\n", "읽을 수 없")):
        with pytest.raises(ValueError, match=why):
            ns.parse_params(bad)
    assert ns.params_sha("filter:\n  min_length: 1\n") == ns.params_sha("filter: {min_length: 1}\n") != ns.params_sha("filter:\n  min_length: 2\n")      # content, not formatting
    ds = datasets.save(1, "Ib", {"clade": ["Ib"]}, include_private=True, build="clade-i")
    ds = datasets.get(ds["id"])
    m = ns.merged_config(ds, "run-1")
    assert m["filter"] == {"min_length": 170000, "min_date": 1900} and m["inputs"][0]["metadata"] == "data/genomecurator/metadata.tsv" and m["use_remote_files"] is False
    assert m["private_metadata"].endswith("metadata-private.tsv") and m["build_name"] == m["auspice_name"] == "run-1"
    datasets.set_params(ds["id"], "filter:\n  min_length: 100\ntimetree: true\n")
    m2 = ns.merged_config(datasets.get(ds["id"]), "run-1")
    assert m2["filter"] == {"min_length": 100, "min_date": 1900} and m2["timetree"] is True and m2["root"] == "best"           # nested keys merge into the workflow's config


@_unix_only
def test_nextstrain_run_keeps_every_run_and_compares_with_today(tmp_path, monkeypatch):
    import json
    import pytest
    from genomecurator import datasets, exclusions, private
    ns = _ns_env(tmp_path, monkeypatch)
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_N1", "strain": "a", "date": "2025", "country": "Uganda", "clade": "Ib"}], {"LAB_N1": "ACGT" * 40}))
    ds = datasets.save(1, "Ib set", {"clade": ["Ib"]}, include_private=True, build="clade-i")
    ds = datasets.get(ds["id"])
    seen = []
    r1 = ns.run(ds, note="first", progress=lambda pct, m: seen.append((round(pct), m)))
    assert r1["status"] == "done" and r1["n_public"] == 1 and r1["n_private"] == 1 and r1["note"] == "first" and r1["build"] == "clade-i" and r1["workflow_commit"] == ""
    assert any("align" in m for _, m in seen) and seen[-1] == (100, "분석 완료") and any("4/4" in m for _, m in seen)
    name = r1["name"]
    assert name.startswith("mpox-ib-set-") and "_" not in name                                    # hyphens only: Auspice reads '_' as a folder separator
    assert (tmp_path / "nextstrain" / "auspice" / f"{name}.json").exists() and (tmp_path / "nextstrain" / "auspice" / f"{name}_root-sequence.json").exists()
    rd = __import__("pathlib").Path(r1["run_dir"])
    assert {p.name for p in rd.iterdir()} >= {"members.txt", "metadata.tsv", "metadata-private.tsv", "manifest.json", "config-overlay.yaml", "run.log", "versions.json"}
    assert ns.run_members(r1) == ["PP_A1", "private:LAB_N1"] and "steps (100%)" in (rd / "run.log").read_text(encoding="utf-8")
    ov = (rd / "config-overlay.yaml").read_text(encoding="utf-8")
    assert f'build_name: "{name}"' in ov and "private_metadata" in ov and "defaults/clade-i/config.yaml" in ov
    # the same members / settings / workflow: not repeated
    again = ns.run(ds)
    assert again["id"] == r1["id"] and again["reused"] and len(ns.list_runs(ds["id"])) == 1
    assert ns.run(ds, force=True)["id"] != r1["id"] and len(ns.list_runs(ds["id"])) == 2          # forced: a new run, the old one stays
    # the data changes -> a new run; the first one is still there and says how it differs from today
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_N2", "strain": "b", "date": "2025", "country": "Uganda", "clade": "Ib"}], {"LAB_N2": "ACGT" * 41}))
    r3 = ns.run(ds, note="after adding LAB_N2")
    assert r3["n_private"] == 2 and not r3.get("reused") and len(ns.list_runs(ds["id"])) == 3
    now = ns.current_members(ds)
    c1 = ns.compare(ns.get_run(r1["id"]), now, ns.params_sha(""), "")
    assert c1["same_members"] is False and c1["added"] == 1 and c1["removed"] == 0 and c1["current"] is False
    c3 = ns.compare(ns.get_run(r3["id"]), now, ns.params_sha(""), "")
    assert c3["current"] and c3["same_members"] and c3["same_params"] and c3["same_workflow"]
    # a settings change alone is a different analysis
    datasets.set_params(ds["id"], "timetree: true\n")
    ds2 = datasets.get(ds["id"])
    assert ns.compare(ns.get_run(r3["id"]), now, ns.params_sha(ds2["params_yaml"]), "")["same_params"] is False
    r4 = ns.run(ds2)
    assert not r4.get("reused") and "timetree: true" in (__import__("pathlib").Path(r4["run_dir"]) / "config-overlay.yaml").read_text(encoding="utf-8")
    # members that are left out change the hash too
    exclusions.add(1, ["PP_A1"], "public", "x")
    assert datasets.members_hash(ns.current_members(ds2)) != r4["members_sha"]
    # a failing workflow is recorded with its log tail, and kept
    monkeypatch.setenv("FAKE_NS_FAIL", "1")
    with pytest.raises(RuntimeError, match="boom"):
        ns.run(ds2, force=True)
    bad = ns.list_runs(ds["id"])[0]
    assert bad["status"] == "error" and "boom" in bad["error"] and not (tmp_path / "nextstrain" / "auspice" / f"{bad['name']}.json").exists()
    monkeypatch.delenv("FAKE_NS_FAIL")
    # an interrupted run (server restart) is marked, notes can be edited, deleting a run removes only its own files
    con = __import__("genomecurator.db", fromlist=["x"]).connect()
    con.execute("UPDATE analysis_run SET status='running' WHERE id=?", (bad["id"],))
    con.commit()
    con.close()
    assert ns.mark_interrupted() == 1 and ns.get_run(bad["id"])["status"] == "interrupted"
    ns.set_note(r1["id"], "keep this")
    assert ns.get_run(r1["id"])["note"] == "keep this"
    assert ns.delete_run(r1["id"]) and not ns.delete_run(r1["id"]) and not (tmp_path / "nextstrain" / "auspice" / f"{r1['name']}.json").exists()
    assert (tmp_path / "nextstrain" / "auspice" / f"{r3['name']}.json").exists() and ns.get_run(r3["id"]) is not None


def test_nextstrain_run_refuses_what_is_missing(tmp_path, monkeypatch):
    import pytest
    from genomecurator import datasets
    ns = _ns_env(tmp_path, monkeypatch)
    ds = datasets.get(datasets.save(1, "Ib", {"clade": ["Ib"]}, build="clade-i")["id"])
    monkeypatch.setattr(ns, "env_prefix", lambda: None)
    with pytest.raises(RuntimeError, match="환경"):
        ns.run(ds)
    monkeypatch.setattr(ns, "env_prefix", lambda: tmp_path / "envs" / "gc-nextstrain")
    monkeypatch.setenv("GC_NEXTSTRAIN_DIR", str(tmp_path / "nowhere"))
    with pytest.raises(RuntimeError, match="워크플로 저장소"):
        ns.run(ds)
    monkeypatch.setenv("GC_NEXTSTRAIN_DIR", str(tmp_path / "wf"))
    with pytest.raises(ValueError, match="빌드"):
        ns.run({**ds, "build": "nope"})
    with pytest.raises(ValueError, match="runtime"):
        ns.run(ds, runtime="cloud")
    with pytest.raises(ValueError, match="멤버가 없습니다"):
        ns.run({**ds, "rules": datasets.normalize({"clade": ["IX"]})})
    assert ns.list_runs(ds["id"]) == []                                                           # nothing was started


@_unix_only
def test_nextstrain_api(tmp_path, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from genomecurator import datasets, web
    ns = _ns_env(tmp_path, monkeypatch)
    ds = datasets.save(1, "Ib set", {"clade": ["Ib"]}, build="clade-i")
    c = TestClient(web.create_app(scheduler=False))
    st = c.get("/api/nextstrain/status").json()
    assert st["env"]["ok"] and st["workflow"]["ok"] and st["workflow"]["builds"] == ["clade-i"] and st["viewer"]["running"] is False
    assert c.get("/nextstrain").status_code == 200 and c.get(f"/nextstrain?dataset={ds['id']}").status_code == 200
    info = c.get(f"/api/datasets/{ds['id']}/nextstrain").json()
    assert info["members"] == 1 and info["qc_missing"] == 0 and "min_length" in info["merged"] and "min_length" in info["defaults"] and info["runs"] == [] and info["error"] == ""
    assert c.get("/api/datasets/999/nextstrain").status_code == 404
    bad = c.post(f"/api/datasets/{ds['id']}/nextstrain/params", json={"params_yaml": "inputs: []"})
    assert bad.status_code == 400 and "관리" in bad.json()["error"]
    assert c.post(f"/api/datasets/{ds['id']}/nextstrain/params", json={"params_yaml": "timetree: true\n"}).json() == {"saved": True}
    assert "timetree: true" in c.get(f"/api/datasets/{ds['id']}/nextstrain").json()["merged"]

    def wait(jid):
        for _ in range(300):
            j = c.get(f"/api/job/{jid}").json()
            if j["done"]:
                return j
            time.sleep(0.05)
        raise AssertionError("job did not finish")
    j = wait(c.post(f"/api/datasets/{ds['id']}/analyze?note=api").json()["job"])
    assert j["status"] == "done" and j["kind"] == "nextstrain" and j["result"]["name"].startswith("mpox-ib-set-") and not j["result"]["reused"]
    j2 = wait(c.post(f"/api/datasets/{ds['id']}/analyze").json()["job"])
    assert j2["result"]["reused"] is True
    runs = c.get(f"/api/analyses?dataset={ds['id']}").json()
    assert len(runs) == 1 and runs[0]["status"] == "done" and runs[0]["url"].endswith("/" + runs[0]["name"]) and runs[0]["note"] == "api"
    d = c.get(f"/api/analyses/{runs[0]['id']}").json()
    assert "steps (100%)" in d["log_tail"] and d["n_public"] == 1
    info = c.get(f"/api/datasets/{ds['id']}/nextstrain").json()
    assert info["runs"][0]["cmp"]["current"] is True
    assert c.post(f"/api/analyses/{runs[0]['id']}/note", json={"note": "good one"}).json() == {"saved": True} and c.get(f"/api/analyses/{runs[0]['id']}").json()["note"] == "good one"
    assert c.get("/api/analyses/999").status_code == 404 and c.post("/api/datasets/999/analyze").status_code == 404
    assert c.delete(f"/api/analyses/{runs[0]['id']}").json() == {"deleted": True} and c.get(f"/api/analyses?dataset={ds['id']}").json() == []
    assert c.post("/api/nextstrain/viewer/restart").status_code == 400 and c.post("/api/nextstrain/viewer/stop").json() == {"stopped": False}


def _sq_env(tmp_path, monkeypatch):
    """_ns_env plus a fake Squirrel environment: `squirrel <fasta> --outdir D ...` writes the files the real one writes and prints snakemake-style progress."""
    import stat
    from genomecurator import squirrel
    ns = _ns_env(tmp_path, monkeypatch)
    env = tmp_path / "envs" / "gc-squirrel"
    (env / "bin").mkdir(parents=True)
    fake = env / "bin" / "squirrel"
    fake.write_text("""#!/usr/bin/env python3
import os, sys
a = sys.argv[1:]
if a[0] in ("-v", "--version"):
    print("squirrel 9.9.9"); sys.exit(0)
fasta, outdir = a[0], a[a.index("--outdir") + 1]
with open(os.path.join(os.path.dirname(os.path.abspath(outdir)), "args.txt"), "a") as f:
    f.write(" ".join(a) + "\\n")
if os.environ.get("FAKE_SQ_FAIL"):
    print("\\x1b[31mnothing aligned\\x1b[0m"); sys.exit(2)
if "-p" in a or "--run-phylo" in a or "--run-apobec3-phylo" in a:
    import subprocess
    subprocess.run(["iqtree3", "-s", "x.fasta"], cwd=os.path.dirname(os.path.abspath(outdir)))     # as Squirrel's own workflow does
names = [l[1:].strip() for l in open(fasta) if l.startswith(">")]
skip = set(l.split(",")[0] for l in open(a[a.index("--exclude") + 1]).read().splitlines()[1:]) if "--exclude" in a else set()
names = [n for n in names if n not in skip]
os.makedirs(outdir, exist_ok=True)
print("rule align:"); print("1 of 3 steps (33%) done"); print("rule tree:"); print("3 of 3 steps (100%) done"); print("Phylogenetics complete.")
open(os.path.join(outdir, "alignment.aln.fasta"), "w").write("".join(f">{n}\\nACGT\\n" for n in names))
open(os.path.join(outdir, "alignment.report.html"), "w").write("<html>report " + ",".join(names) + "</html>")
if ("-p" in a or "--run-phylo" in a or "--run-apobec3-phylo" in a) and not os.environ.get("FAKE_SQ_NOTREE"):
    open(os.path.join(outdir, "alignment.tree"), "w").write("#NEXUS\\n")
open(os.path.join(outdir, "alignment.suggested_mask.csv"), "w").write("Name,Minimum\\n")
ex = os.environ.get("FAKE_SQ_EXCLUDE")
open(os.path.join(outdir, "suggested_to_exclude.csv"), "w").write("name,note\\n" + (f"{ex},N content is 0.3\\n" if ex else ""))
""")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    it = env / "bin" / "iqtree"
    it.write_text("#!/bin/sh\necho 'IQ-TREE version 9.9 built today'\n")
    it.chmod(it.stat().st_mode | stat.S_IEXEC)
    it3 = env / "bin" / "iqtree3"
    it3.write_text("#!/bin/sh\necho \"$@\" >> iq_args.txt\n")
    it3.chmod(it3.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setattr(squirrel, "env_prefix", lambda: env)
    return squirrel


def test_squirrel_settings_and_command(tmp_path, monkeypatch):
    import pytest
    from genomecurator import datasets
    sq = _sq_env(tmp_path, monkeypatch)
    ds = datasets.get(datasets.save(1, "Ib", {"clade": ["Ib"]}, build="clade-i")["id"])
    assert sq.defaults_for(ds)["clade"] == "cladei" and sq.defaults_for({"build": "hmpxv1"})["clade"] == "cladeiib" and sq.defaults_for({})["clade"] == "split"
    assert sq.parse_settings("") == {} and sq.parse_settings({"clade": "cladeib", "seq_qc": 0, "threads": "8", "extra_args": "--no-mask --outgroups KJ642615"})["threads"] == 8
    for bad, why in (({"clade": "cladeX"}, "clade"), ({"bogus": 1}, "알 수 없는"), ({"threads": "0"}, "threads"), ({"extra_args": "--outdir /etc"}, "허용되지"),
                     ({"extra_args": "--tempdir x"}, "허용되지"), ({"extra_args": "--outgroups ../x"}, "허용되지"), ({"extra_args": "--no-mask; rm -rf /"}, "허용되지")):
        with pytest.raises(ValueError, match=why):
            sq.parse_settings(bad)
    with pytest.raises(ValueError, match="JSON"):
        sq.parse_settings("{not json")
    cmd = sq.command(ds, "in.fa", "out", "tmp", 4)
    assert cmd[:2] == ["squirrel", "in.fa"] and cmd[cmd.index("--clade") + 1] == "cladei" and "--seq-qc" in cmd and "--run-phylo" in cmd and "--include-background" not in cmd and cmd[cmd.index("-t") + 1] == "4"
    ds2 = {**ds, "squirrel_json": '{"apobec3": true, "include_background": true, "seq_qc": false, "extra_args": "--no-mask"}'}
    c2 = sq.command(ds2, "in.fa", "out", "tmp", 2)
    assert "--run-apobec3-phylo" in c2 and "--run-phylo" not in c2 and "--include-background" in c2 and "--seq-qc" not in c2 and c2[-1] == "--no-mask"
    assert sq.settings_sha(ds) == sq.settings_sha({**ds, "squirrel_json": '{"clade": "cladei"}'}) != sq.settings_sha(ds2)       # the effective settings, not how they were written


@_unix_only
def test_squirrel_run_keeps_runs_beside_nextstrain_ones(tmp_path, monkeypatch):
    import json
    import pytest
    from pathlib import Path
    from genomecurator import datasets, private
    sq = _sq_env(tmp_path, monkeypatch)
    ns = __import__("genomecurator.nextstrain", fromlist=["x"])
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_S1", "strain": "a", "date": "2025", "country": "Uganda", "clade": "Ib"}], {"LAB_S1": "ACGT" * 40}))
    ds = datasets.get(datasets.save(1, "Ib set", {"clade": ["Ib"]}, include_private=True, build="clade-i")["id"])
    seen = []
    r1 = sq.run(ds, note="first", progress=lambda pct, m: seen.append((round(pct), m)))
    assert r1["status"] == "done" and r1["tool"] == "squirrel" and r1["n_public"] == 1 and r1["n_private"] == 1 and "-squirrel-" in r1["name"] and r1["workflow_commit"] == "squirrel 9.9.9"
    assert any("align" in m for _, m in seen) and seen[-1] == (100, "분석 완료")
    rd = Path(r1["run_dir"])
    args = (rd / "args.txt").read_text(encoding="utf-8")
    passes = args.strip().splitlines()
    assert len(passes) == 2 and "--run-phylo" not in passes[0] and "--seq-qc" in passes[0] and "--run-phylo" in passes[1] and "--clade cladei" in passes[1]     # QC first, then the tree
    assert "--exclude" not in passes[1]                                                      # nothing was flagged: nothing to leave out
    assert (rd / "iq_args.txt").read_text().strip() == "-st DNA -s x.fasta"                  # IQ-TREE was told the data type (and only in the tree pass)
    assert not (rd / "shim").exists() and "-st DNA" in (rd / "run.log").read_text(encoding="utf-8")
    assert r1["auspice_file"] == "squirrel/alignment.report.html"
    rep = (rd / "squirrel" / "alignment.report.html").read_text(encoding="utf-8")
    assert "PP_A1" in rep and "LAB_S1" in rep                                               # the dataset's public AND private members went in
    assert (rd / "squirrel" / "alignment.aln.fasta.gz").exists() and not (rd / "squirrel" / "alignment.aln.fasta").exists() and not (rd / "input.fasta").exists() and not (rd / "tmp").exists()
    assert ns.run_members(r1) == ["PP_A1", "private:LAB_S1"] and "steps (100%)" in (rd / "run.log").read_text(encoding="utf-8")
    assert {f["name"] for f in sq.result_files(r1)} >= {"alignment.report.html", "alignment.tree", "alignment.suggested_mask.csv", "alignment.aln.fasta.gz"}
    again = sq.run(ds)
    assert again["id"] == r1["id"] and again["reused"]
    # a Nextstrain run of the same dataset is a different thing: it neither blocks nor is blocked
    n1 = ns.run(ds)
    assert n1["tool"] == "nextstrain" and not n1.get("reused") and ns.run(ds)["reused"] and sq.run(ds)["reused"]
    assert [r["tool"] for r in ns.list_runs(ds["id"])] == ["nextstrain", "squirrel"]
    # settings changed -> a new analysis; forced -> too
    datasets.set_squirrel(ds["id"], '{"apobec3": true}')
    ds2 = datasets.get(ds["id"])
    r2 = sq.run(ds2)
    assert not r2.get("reused") and r2["id"] != r1["id"] and "--run-apobec3-phylo" in (Path(r2["run_dir"]) / "args.txt").read_text()
    # what Squirrel's own QC flags (N content ...) is left out of the tree and reported
    monkeypatch.setenv("FAKE_SQ_EXCLUDE", "PP_A1")
    r5 = sq.run(ds2, force=True)
    p5 = (Path(r5["run_dir"]) / "args.txt").read_text().strip().splitlines()
    assert "--exclude" in p5[1] and "PP_A1" not in (Path(r5["run_dir"]) / "squirrel" / "alignment.report.html").read_text()
    assert (Path(r5["run_dir"]) / "squirrel" / "qc_suggested_to_exclude.csv").exists() and json.loads((Path(r5["run_dir"]) / "versions.json").read_text())["excluded_by_qc"] == ["PP_A1"]
    assert "제외 1건" in (Path(r5["run_dir"]) / "run.log").read_text() and not (Path(r5["run_dir"]) / "squirrel_qc").exists()
    datasets.set_squirrel(ds["id"], '{"apobec3": true, "exclude_flagged": false}')
    r6 = sq.run(datasets.get(ds["id"]))
    assert "--exclude" not in (Path(r6["run_dir"]) / "args.txt").read_text().strip().splitlines()[1]                       # the person can keep everything
    monkeypatch.delenv("FAKE_SQ_EXCLUDE")
    # without a phylogeny it is a single call
    datasets.set_squirrel(ds["id"], '{"run_phylo": false, "apobec3": false}')
    r7 = sq.run(datasets.get(ds["id"]))
    assert len((Path(r7["run_dir"]) / "args.txt").read_text().strip().splitlines()) == 1 and not list((Path(r7["run_dir"]) / "squirrel").glob("*.tree"))
    # asked for a tree but Squirrel exited 0 without one (what IQ-TREE did on the real data): an error, not a "done"
    datasets.set_squirrel(ds["id"], '{"apobec3": true}')
    monkeypatch.setenv("FAKE_SQ_NOTREE", "1")
    with pytest.raises(RuntimeError, match="계통수를 만들지 못했습니다"):
        sq.run(datasets.get(ds["id"]), force=True)
    assert ns.list_runs(ds["id"])[0]["status"] == "error"
    monkeypatch.delenv("FAKE_SQ_NOTREE")
    datasets.set_squirrel(ds["id"], '{"apobec3": true}')
    ds2 = datasets.get(ds["id"])
    assert sq.run(ds2, force=True)["id"] != r2["id"]
    # failure: recorded with the (colour-stripped) log tail
    monkeypatch.setenv("FAKE_SQ_FAIL", "1")
    with pytest.raises(RuntimeError, match="nothing aligned"):
        sq.run(ds2, force=True)
    bad = ns.list_runs(ds["id"])[0]
    assert bad["tool"] == "squirrel" and bad["status"] == "error" and "\x1b" not in bad["error"]


def test_squirrel_refuses_what_is_missing(tmp_path, monkeypatch):
    import pytest
    from genomecurator import datasets
    sq = _sq_env(tmp_path, monkeypatch)
    ds = datasets.get(datasets.save(1, "Ib", {"clade": ["Ib"]}, build="clade-i")["id"])
    monkeypatch.setattr(sq, "env_prefix", lambda: None)
    with pytest.raises(RuntimeError, match="환경"):
        sq.run(ds)
    monkeypatch.setattr(sq, "env_prefix", lambda: tmp_path / "envs" / "gc-squirrel")
    with pytest.raises(ValueError, match="멤버가 없습니다"):
        sq.run({**ds, "rules": datasets.normalize({"clade": ["IX"]})})


@_unix_only
def test_squirrel_api_and_files(tmp_path, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from genomecurator import datasets, web
    sq = _sq_env(tmp_path, monkeypatch)
    ds = datasets.save(1, "Ib set", {"clade": ["Ib"]}, build="clade-i")
    c = TestClient(web.create_app(scheduler=False))
    st = c.get("/api/nextstrain/status").json()
    assert st["squirrel"]["ok"] and st["squirrel"]["version"] == "squirrel 9.9.9"
    info = c.get(f"/api/datasets/{ds['id']}/nextstrain").json()
    assert info["squirrel"]["settings"]["clade"] == "cladei" and "split" in info["squirrel"]["clades"]
    bad = c.post(f"/api/datasets/{ds['id']}/squirrel/settings", json={"settings": {"extra_args": "--outdir /x"}})
    assert bad.status_code == 400 and c.post("/api/datasets/999/squirrel/settings", json={"settings": {}}).status_code == 404
    assert c.post(f"/api/datasets/{ds['id']}/squirrel/settings", json={"settings": {"clade": "cladeib", "apobec3": True}}).json()["settings"] == {"clade": "cladeib", "apobec3": True}
    assert c.get(f"/api/datasets/{ds['id']}/nextstrain").json()["squirrel"]["settings"]["clade"] == "cladeib"

    def wait(jid):
        for _ in range(300):
            j = c.get(f"/api/job/{jid}").json()
            if j["done"]:
                return j
            time.sleep(0.05)
        raise AssertionError("job did not finish")
    j = wait(c.post(f"/api/datasets/{ds['id']}/squirrel?note=api").json()["job"])
    assert j["status"] == "done" and j["kind"] == "squirrel" and not j["result"]["reused"]
    assert wait(c.post(f"/api/datasets/{ds['id']}/squirrel").json()["job"])["result"]["reused"] is True
    runs = c.get(f"/api/datasets/{ds['id']}/nextstrain").json()["runs"]
    r = runs[0]
    assert r["tool"] == "squirrel" and r["note"] == "api" and r["cmp"]["current"] is True and r["url"].endswith("alignment.report.html") and {f["name"] for f in r["files"]} >= {"alignment.tree"}
    rep = c.get(r["url"])
    assert rep.status_code == 200 and "report PP_A1" in rep.text
    assert c.get(f"/api/analyses/{r['id']}/files/squirrel/alignment.tree").status_code == 200
    assert c.get(f"/api/analyses/{r['id']}/files/..%2F..%2Fgenomecurator.db").status_code == 404 and c.get(f"/api/analyses/{r['id']}/files/squirrel/nope.txt").status_code == 404
    assert c.get(f"/api/analyses/{r['id']}/files/../../g.db").status_code in (404, 400, 422) and c.get("/api/analyses/999/files/x").status_code == 404
    assert c.post("/api/datasets/999/squirrel").status_code == 404


def test_mutprofile_annotation_and_paf():
    from genomecurator import mutprofile as mp
    # a plus-strand gene 4..12 (ATG AAA TGA) and a minus-strand gene 20..28 whose coding sequence is ATG GGG TAA
    ref = "CCCATGAAATGACCCCCCCCCCATGGGGTAACC".replace("CCCCCCCCCC", "CCCCCCC")
    ref = ref[:3] + "ATGAAATGA" + "C" * 7 + mp_rev("ATGGGGTAA") + "CC"
    cds = [{"gene": "G1", "product": "p1", "strand": "+", "exons": [(4, 12)]}, {"gene": "G2", "product": "p2", "strand": "-", "exons": [(20, 28)]}]
    a = mp.annotate(8, "A", "G", cds, ref)                         # AAA -> AGA: K2R, nonsynonymous
    assert (a["gene"], a["aa_change"], a["effect"], a["codon_base"]) == ("G1", "K2R", "nonsynonymous", 2)
    assert mp.annotate(9, "A", "G", cds, ref)["effect"] == "synonymous"   # AAA -> AAG: K -> K
    assert mp.annotate(1, "C", "T", cds, ref)["effect"] == "intergenic"
    b = mp.annotate(27, ref[26], "T", cds, ref)                    # minus strand: plus base A (coding T of ATG) -> T (coding A): ATG -> AAG
    assert (b["gene"], b["strand"], b["aa_change"], b["codon_base"]) == ("G2", "-", "M1K", 2)
    paf = ["q1\t100\t0\t100\t+\tref\t1000\t10\t110\t97\t100\t60\ttp:A:P\tcs:Z::5*ag:3*nt:2+ac:4-gg:10*ct:20"]
    assert list(mp.snvs_from_paf(paf)) == [("q1", 16, "A", "G"), ("q1", 39, "C", "T")]      # the N substitution (*nt) is never reported; insertions/deletions move only the right coordinate


def test_mutprofile_reference_choice():
    from genomecurator import mutprofile as mp
    assert mp.reference_choice({"build": "clade-i"}) == "KJ642613.1"                      # clade I follows the Nextstrain clade-i build's reference
    assert mp.reference_choice({"build": "mpxv"}) == mp.reference_choice({"build": ""}) == "NC_063383.1"
    assert mp.reference_choice({"build": "clade-i", "mutprofile_json": '{"reference": "NC_063383.1"}'}) == "NC_063383.1"
    assert mp.reference_choice({"build": "clade-i", "mutprofile_json": '{"reference": "NC_003310.1"}'}) == "NC_003310.1"
    assert mp.settings_sha({"build": "clade-i"}) != mp.settings_sha({"build": "mpxv"})   # a different reference is a different analysis


def test_genbank_reference_to_fasta_and_gff3():
    from genomecurator import mutprofile as mp
    gb = '''LOCUS       TEST       30 bp    DNA     linear   VRL 01-JAN-2000
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="x"
     CDS             join(1..6,
                     10..15)
                     /gene="A1"
                     /product="two-exon;
                     protein"
     CDS             complement(20..28)
                     /gene="B2"
ORIGIN
        1 atgaaatttc atgcccccca tgggggtaac
//
'''
    seq, gff = mp.genbank_to_fasta_gff3(gb, "TEST")
    assert len(seq) == 30 and seq.startswith("ATGAAA")
    rows = [l.split("\t") for l in gff.splitlines() if not l.startswith("#")]
    assert [(r[3], r[4], r[6]) for r in rows] == [("1", "6", "+"), ("10", "15", "+"), ("20", "28", "-")]
    assert "gene=A1" in rows[0][8] and "product=two-exon%3B protein" in rows[0][8] and rows[0][8].split(";")[0] == rows[1][8].split(";")[0]


def mp_rev(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_pathoplexus_prepare_package_never_uploads(tmp_path, monkeypatch):
    import json, zipfile, io
    from pathlib import Path
    from genomecurator import pathoplexus, private
    _qc_env(tmp_path, monkeypatch)
    monkeypatch.setattr(pathoplexus, "template", lambda org, refresh=False: list(pathoplexus.FALLBACK_TEMPLATE))
    import urllib.request
    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **k: (_ for _ in ()).throw(AssertionError("the app must not touch the network here")))
    cols = ["PPX_accession", "strain", "date", "country", "division", "host", "institution", "authors"]
    rws = [["LAB_001", "s1", "2024-03-05", "Uganda", "Kampala", "Homo sapiens", "Makerere University", "Kim, Min; Lee, Jo"], ["LAB_002", "s2", "2024", "Uganda", "", "Homo sapiens", "Makerere University", "no comma here"]]
    files = ("\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rws), ">LAB_001\n" + "ACGT" * 50 + "\n>LAB_002\n" + "ACGT" * 40 + "\n")
    private.import_texts(1, *files, source_name="b.tsv")
    c = pathoplexus.candidates(1)
    rows = {r["id"]: r for r in c["rows"]}
    assert rows["LAB_001"]["status"] in ("ok", "warn") and not rows["LAB_001"]["errors"]
    assert any('authors' in m for m in rows["LAB_002"]["errors"])           # authors must be 'Last, First; ...'
    with __import__("pytest").raises(ValueError):
        pathoplexus.build_package(1, ["LAB_001"], "RESTRICTED", "")  # a restricted package needs an end date
    p = pathoplexus.build_package(1, ["LAB_001"], "OPEN", "", "42", "test")
    d = Path(p["dir"])
    assert {f.name for f in d.iterdir()} >= {"metadata.tsv", "sequences.fasta", "submission.json", "submit.sh", "README.txt"}
    meta = (d / "metadata.tsv").read_text(encoding="utf-8").splitlines()
    assert meta[0].split("\t")[0] == "id" and meta[1].split("\t")[0] == "LAB_001" and len(meta) == 2
    assert (d / "sequences.fasta").read_text(encoding="utf-8").startswith(">LAB_001")
    sh = (d / "submit.sh").read_text(encoding="utf-8")
    assert "PATHOPLEXUS_TOKEN" in sh and "eyJ" not in sh                # the token is read from the environment, never stored
    assert pathoplexus.mark(p["id"], "submitted")["status"] == "submitted"
    assert zipfile.ZipFile(io.BytesIO(pathoplexus.zip_package(p["id"]))).namelist()
    assert {r["id"]: r for r in pathoplexus.candidates(1)["rows"]}["LAB_001"]["submitted"]["status"] == "submitted"
    # the page and API
    from fastapi.testclient import TestClient
    from genomecurator import web
    cl = TestClient(web.create_app(scheduler=False))
    assert cl.get("/pathoplexus").status_code == 200 and cl.get("/api/pathoplexus/1/candidates").json()["rows"]
    assert cl.get("/api/pathoplexus/packages").json()[0]["id"] == p["id"]


def test_snipit_data_filters_groups_and_hides_common(tmp_path, monkeypatch):
    from genomecurator import mutprofile as mp
    rd = tmp_path / "run"
    (rd / "mutprofile").mkdir(parents=True)
    hdr = "sample\tpos\tref\talt\tgene\tproduct\tstrand\taa_change\teffect\tref_codon\talt_codon\tcodon_pos\tcodon_base\tmotif_ref\tobserved_profile\tprofile_of_interest\n"
    def row(s, pos, alt, poi): return f"{s}\t{pos}\tG\t{alt}\tG1\tp\t+\tK2R\tnonsynonymous\t\t\t\t\t\t\t{poi}\n"
    (rd / "mutprofile" / "mutations.tsv").write_text(hdr + row("a", 10, "A", "GA>AA") + row("a", 50, "T", "other") + row("b", 10, "A", "GA>AA") + row("c", 10, "A", "GA>AA") + row("c", 60, "C", "other") + row("c", 70, "C", "other"), encoding="utf-8")
    (rd / "mutprofile" / "summary.json").write_text('{"reference": "R"}', encoding="utf-8")
    ref = tmp_path / "r.fasta"
    ref.write_text(">R\n" + "ACGT" * 50 + "\n", encoding="utf-8")
    monkeypatch.setattr(mp, "reference_by_name", lambda acc: (ref, None))
    monkeypatch.setattr(mp, "_run_meta", lambda run: {"a": {"clade": "Ib", "lineage": "A.1", "country": "", "date": ""}, "b": {"clade": "Ib", "lineage": "A.1", "country": "", "date": ""}, "c": {"clade": "Ib", "lineage": "A.2", "country": "", "date": ""}})
    run = {"run_dir": str(rd)}
    d = mp.snipit_data(run)
    assert d["reference"] == "R" and d["length"] == 200 and [s["id"] for s in d["samples"]] == ["a", "b", "c"] and d["lineages"] == ["A.1", "A.2"]
    assert [s["id"] for s in mp.snipit_data(run, per_group=1)["samples"]] == ["a", "c"]               # one per lineage, the one with most SNPs first
    assert [s["id"] for s in mp.snipit_data(run, lineages=["A.2"])["samples"]] == ["c"]
    assert all(m[6] for s in mp.snipit_data(run, apobec_only=True)["samples"] for m in s["snps"])
    assert all(m[0] != 10 for s in mp.snipit_data(run, hide_common=1.0)["samples"] for m in s["snps"])  # present in every shown sample -> hidden
    assert all(55 <= m[0] <= 65 for s in mp.snipit_data(run, start=55, end=65)["samples"] for m in s["snps"])


def test_lineage_definitions_chain_dates_and_spread(tmp_path, monkeypatch):
    from genomecurator import lineages
    base = tmp_path / "ld" / "definitions" / "Ib"
    (base / "sh2023").mkdir(parents=True)
    (base / "sh2023.yml").write_text("name: sh2023\ndefining_snps:\n  - pos: 100\n    nucleotide: C\n", encoding="utf-8")
    (base / "sh2023" / "A.yml").write_text("name: A\nparent: ''\ndesignation_date: '2024-01-01'\ndefining_snps:\n  - pos: 100\n    nucleotide: C\n  - pos: 200\n    nucleotide: G\n", encoding="utf-8")
    (base / "sh2023" / "A.1.yml").write_text("name: A.1\nparent: A\ndesignation_date: '2024-06-01'\ndefining_snps:\n  - pos: 300\n    nucleotide: T\n", encoding="utf-8")
    monkeypatch.setenv("GC_LINEAGE_DIR", str(tmp_path / "ld"))
    lineages._DEF.clear()
    assert [(p, n, by) for p, n, by in lineages.defining_snps("sh2023", "A.1")] == [(100, "C", "sh2023"), (200, "G", "A"), (300, "T", "A.1")]   # stepwise, a repeated SNP kept once
    assert lineages.defining_snps("sh2023", "nope") == [(100, "C", "sh2023")] and lineages.status()["ok"]
    assert lineages.parse_date("2024-11-21")[:3] == ("day", "2024-W47", "2024-11") and lineages.parse_date("2024-12-30")[1] == "2025-W01"      # ISO year differs from the calendar year
    assert lineages.parse_date("2025-01-XX")[:3] == ("month", None, "2025-01") and lineages.parse_date("2025-XX-XX")[0] == "" and lineages.parse_date("2025-02-31")[0] == ""
    gs = [{"id": f"g{i}", "lineage": lin, "outbreak": "sh2023", "clade": "Ib", "country": c, "date": d} for i, (lin, c, d) in enumerate(
        [("A", "Uganda", "2024-02-05"), ("A", "Uganda", "2024-02-06"), ("A", "Kenya", "2024-02-20"), ("A", "Kenya", "2024-02-21"), ("A", "Kenya", "2024-03-XX"), ("A.1", "Uganda", "2024-07-01"), ("A.1", "Uganda", "2024-07-02"), ("A.1", "Uganda", "2024-07-03"), ("A.1", "Germany", "2024-08-XX"), ("A.1", "Germany", "")])]
    monkeypatch.setattr(lineages, "_genomes", lambda run: [dict(g) for g in gs])
    monkeypatch.setattr(lineages.mutprofile, "summary", lambda run: {"reference": "R"})
    d = lineages.spread({}, min_n=1)
    assert d["bin"] == "week" and d["no_date"] == 1 and d["month_only"] == 2 and [r["label"] for r in d["rows"]] == ["A", "A.1"]
    a = d["rows"][0]
    assert (a["first"], a["last"], a["n"]) == ("2024-02-05", "2024-03", 5) and a["own"] == ["100C", "200G"] and a["snps"] == ["100C", "200G"]
    kenya = next(c for c in a["countries"] if c["country"] == "Kenya")
    assert kenya["cells"] == {"2024-W08": 2} and kenya["monthly"] == {"2024-03": 1}                    # a month-only date is a faint month band, never a week
    assert lineages.spread({}, bin_="month", min_n=1)["bin"] == "month" and [r["label"] for r in lineages.spread({}, min_n=5)["rows"]] == ["A"] and lineages.spread({}, min_n=9)["rows"] == []
    sparse = [{"id": f"s{i}", "lineage": "A", "outbreak": "sh2023", "clade": "", "country": "U", "date": f"2024-{(i % 12) + 1:02d}-{(i % 3) * 9 + 1:02d}"} for i in range(12)]
    monkeypatch.setattr(lineages, "_genomes", lambda run: sparse)
    assert lineages.spread({}, min_n=1)["bin"] == "month"                                                # under two genomes per occupied week -> month


def test_insight_apobec_sites_branches_emerging_and_signatures(tmp_path, monkeypatch):
    import json
    from genomecurator import insight, lineages, mutprofile, nextstrain
    t = insight.chi2_2x2(50, 50, 10, 90)
    assert t["chi2"] > 35 and t["p"] < 0.001 and t["odds_ratio"] == 9.0 and insight.chi2_2x2(0, 5, 0, 5)["chi2"] is None
    assert insight._apo_branch("G", "A", 3, "CCGAT") and not insight._apo_branch("G", "A", 3, "CCGTT") and insight._apo_branch("C", "T", 4, "CCTCT") and not insight._apo_branch("C", "T", 4, "CCACT")
    hdr = "sample\tpos\tref\talt\tgene\tproduct\tstrand\taa_change\teffect\tref_codon\talt_codon\tcodon_pos\tcodon_base\tmotif_ref\tobserved_profile\tprofile_of_interest\n"
    rows = []
    def add(s, pos, alt, apo): rows.append(f"{s}\t{pos}\tG\t{alt}\tG1\tp\t+\tK2R\tnonsynonymous\t\t\t\t\t\t\t{'GA>AA' if apo else 'other'}\n")
    for s in "abcd":
        add(s, 10, "A", True)                                    # carried by all four -> fixed (> 50 % of genomes would be 'fixed' only above 50 carriers; here a shared site)
    for s in "ab":
        add(s, 20, "T", False)
    add("a", 30, "A", True); add("b", 31, "C", False); add("c", 32, "T", True)
    rd = tmp_path / "run"
    (rd / "mutprofile").mkdir(parents=True)
    (rd / "mutprofile" / "mutations.tsv").write_text(hdr + "".join(rows), encoding="utf-8")
    run = {"run_dir": str(rd), "dataset_id": 1}
    monkeypatch.setattr(mutprofile, "summary", lambda r: {"reference": ""})
    a = insight.apobec_by_sites(run)
    assert a["private"]["sites"] == 3 and a["private"]["apobec_sites"] == 2 and a["shared"]["sites"] == 2 and a["shared"]["apobec_sites"] == 1
    assert [c["carriers"] for c in a["classes"]] == ["1", "2-3", "4-10", "11-50", ">50"] and a["classes"][1]["sites"] == 1 and a["classes"][2]["sites"] == 1
    monkeypatch.setattr(nextstrain, "list_runs", lambda *a, **k: [])
    assert insight.apobec_by_branches(run) is None                                                       # no tree -> no branch test
    # emerging: a mutation carried by 3 recent genomes, not a defining one, is listed and 'recent'
    gs = [{"id": s, "lineage": "A", "outbreak": "sh2023", "clade": "Ib", "country": c, "date": d} for s, c, d in [("a", "X", "2026-01-01"), ("b", "Y", "2026-02-01"), ("c", "Y", "2026-03-01"), ("d", "Z", "2020-01-01")]]
    monkeypatch.setattr(lineages, "_genomes", lambda r: gs)
    monkeypatch.setattr(lineages, "status", lambda: {"ok": False, "commit": ""})
    e = insight.emerging(run, window_days=100, min_carriers=2, max_share=1.0)
    byk = {(m["pos"], m["alt"]): m for m in e["mutations"]}
    assert (20, "T") in byk and byk[(20, "T")]["n"] == 2 and byk[(20, "T")]["countries"] == ["X", "Y"] and (30, "A") not in byk       # a single carrier is not listed
    assert byk[(20, "T")]["recent"] and e["definitions"] is False
    assert any(set(c["mutations"]) >= {"G10A", "G20T"} for c in e["clusters"]) is False                    # different carrier sets do not form a cluster


def test_outbreak_statistics_clusters_and_categories(monkeypatch):
    from genomecurator import lineages, outbreak
    d = outbreak.size_distribution([53, 42, 7, 10, 6, 11, 3])                 # sizes like the paper's seven domestic subclusters
    assert abs(d["vmr"] - 21.2) < 0.3 and d["vmr"] == d["pearson"] and d["nb_better"] and d["overdispersed"]
    assert d["aic_nb"] < d["aic_poisson"] - 50 and 1.0 < d["theta"] < 1.6 and outbreak.size_distribution([5, 5, 5, 5])["theta"] > 100 and outbreak.size_distribution([1, 2]) is None
    t = outbreak.trend([2022, 2023, 2024, 2025, 2026], [10, 20, 30, 41, 50])
    assert abs(t["slope_per_year"] - 10.1) < 0.01 and t["rho"] == 1.0
    cat = outbreak.categories()
    assert cat["OPG056"]["category"] == "Assembly, budding" and cat["OPG153"]["category"] == "Surface protein" and cat["OPG001"]["category"] == "Host-modulation"
    assert outbreak._sites("ATG") == 0 and abs(outbreak._sites("TTT") - 1 / 3) < 1e-9 and outbreak._sites("TAA") == 0
    # a small tree: root -> X (G10A; 6 tips: Y has 5 of them with a second branch mutation T20C) and Z (2 tips)
    def tip(n): return {"name": n}
    ytips = [tip(f"y{i}") for i in range(5)]
    tree = {"name": "root", "children": [
        {"name": "X", "branch_attrs": {"mutations": {"nuc": ["G10A"]}}, "children": [{"name": "Y", "branch_attrs": {"mutations": {"nuc": ["T20C", "G30A"]}}, "children": ytips}, tip("x1")]},
        {"name": "Z", "branch_attrs": {"mutations": {"nuc": ["G10A"]}}, "children": [tip("z1"), tip("z2")]}]}
    ref = "A" * 9 + "G" + "A" * 9 + "T" + "A" * 9 + "G" + "A" * 20
    monkeypatch.setattr(outbreak, "_tree", lambda run: (tree, ref, {"name": "t", "build": ""}))
    monkeypatch.setattr(outbreak, "_tree_reference", lambda ns: ("R", ref, []))
    gs = [{"id": i, "lineage": "A", "outbreak": "o", "clade": "Ib", "country": "U" if i.startswith("y") else "K", "date": f"2024-0{1 + n % 3}-05"} for n, i in enumerate(["y0", "y1", "y2", "y3", "y4", "x1", "z1", "z2"])]
    monkeypatch.setattr(lineages, "_genomes", lambda run: gs)
    c = outbreak.clusters({"dataset_id": 1}, min_size=5)
    assert [x["n"] for x in c["clusters"]] == [5] and c["unassigned"] == 3 and c["clusters"][0]["countries"] == {"U": 5}
    labs = {m["label"]: m for m in c["clusters"][0]["mutations"]}
    assert set(labs) == {"T20C", "G30A"} and c["clusters"][0]["first"] == "2024-01-05"
    assert not labs["T20C"]["homoplasy"] and [x["label"] for x in c["columns"]] == ["T20C", "G30A"]
    assert c["matrix"][0]["cells"] == {"T20C": "C", "G30A": "A"} and "G10A" in c["clusters"][0]["path"]   # the parent branch's mutation is part of the cluster's path


def test_cross_site_writes_are_refused_for_every_route(tmp_path, monkeypatch):
    _qc_env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    c = TestClient(web.create_app(scheduler=False))
    for path in ("/deploy", "/changeset/1/apply", "/watch/x/check", "/genomes/fasta", "/api/datasets/preview"):      # routes that had no check of their own
        r = c.post(path, headers={"origin": "https://evil.example"})
        assert r.status_code == 403 and r.json()["error"] == "cross-site request refused", path
        assert c.post(path, headers={"origin": "null"}).status_code == 403
        assert c.post(path, headers={"sec-fetch-site": "cross-site"}).status_code == 403
    r = c.post("/api/datasets/preview", json={}, headers={"origin": "http://testserver"})                         # same origin is not stopped by the guard
    assert r.status_code != 403
    assert c.get("/", headers={"origin": "https://evil.example"}).status_code == 200                              # reading stays possible


@_unix_only
def test_changed_content_under_the_same_ids_is_not_the_same_analysis(tmp_path, monkeypatch):
    from genomecurator import datasets, nextstrain, private
    from genomecurator.db import connect
    sq = _sq_env(tmp_path, monkeypatch)
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_C1", "strain": "a", "date": "2025", "country": "Uganda", "clade": "Ib"}], {"LAB_C1": "ACGT" * 40}))
    ds = datasets.get(datasets.save(1, "Ib set", {"clade": ["Ib"]}, include_private=True, build="clade-i")["id"])
    r1 = sq.run(ds)
    assert sq.run(ds)["id"] == r1["id"]                                                      # unchanged -> reused
    base = nextstrain.content_sha(ds)
    cmp0 = nextstrain.compare(nextstrain.get_run(r1["id"]), nextstrain.current_members(ds), sq.settings_sha(ds), r1["workflow_commit"], base)
    assert cmp0["same_content"] is True and cmp0["current"] is True
    con = connect()                                                                           # a curated value of a public member changes (same ids)
    con.execute("UPDATE genome SET isolation_country='Kenya', clade='Ib' WHERE accession='A1'")
    con.commit()
    con.close()
    changed = nextstrain.content_sha(ds)
    r_ = nextstrain.get_run(r1["id"])
    cmp1 = nextstrain.compare(r_, nextstrain.current_members(ds), sq.settings_sha(ds), r1["workflow_commit"], changed)
    assert changed != base
    assert cmp1["same_content"] is False and cmp1["current"] is False
    assert sq.run(ds)["id"] != r1["id"]                                                      # and a new run is made, not the old one reused
    # a private genome gets a new version -> different content
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_C1", "strain": "a", "date": "2025", "country": "Kenya", "clade": "Ib"}], {"LAB_C1": "ACGT" * 40}), source_name="v2")
    assert nextstrain.content_sha(ds) != changed
    old = dict(r_); old["content_sha"] = None                                                 # a run from before the fingerprint existed
    assert nextstrain.compare(old, nextstrain.current_members(ds), "x", "y", changed)["same_content"] is None


def test_itr_mask_unknown_bases_and_within_lineage_trend(tmp_path, monkeypatch):
    import json
    from genomecurator import lineages, mutprofile, outbreak
    assert mutprofile.in_itr(100, 196442) and mutprofile.in_itr(196000, 196442) and not mutprofile.in_itr(100000, 196442) and not mutprofile.in_itr(10, 100)   # a tiny reference is never masked
    assert mutprofile.parse_settings({"mask_itr": False, "itr_bp": "3000"}) == {"mask_itr": False, "itr_bp": 3000}
    with __import__("pytest").raises(ValueError):
        mutprofile.parse_settings({"itr_bp": 99999})
    # ITR variants are left out when the run's mutations are read (also for runs made before the setting existed)
    rd = tmp_path / "run" / "mutprofile"
    rd.mkdir(parents=True)
    hdr = "sample\tpos\tref\talt\tgene\tproduct\tstrand\taa_change\teffect\tref_codon\talt_codon\tcodon_pos\tcodon_base\tmotif_ref\tobserved_profile\tprofile_of_interest\n"
    row = lambda s, p: f"{s}\t{p}\tG\tA\t\t\t\t\tintergenic\t\t\t\t\t\t\tother\n"
    (rd / "mutations.tsv").write_text(hdr + row("a", 100) + row("a", 50000) + row("a", 196000), encoding="utf-8")
    (rd / "summary.json").write_text(json.dumps({"reference": "R", "reference_length": 196442}), encoding="utf-8")
    run = {"run_dir": str(tmp_path / "run"), "dataset_id": 1}
    assert [int(r["pos"]) for r in mutprofile._load_mutations(run)] == [50000]
    assert mutprofile.run_itr(run) == (6501, 189942)
    (rd / "samples.tsv").write_text("id\tlength\tn_frac\na\t196000\t0.01\nb\t150000\t0.4\n", encoding="utf-8")
    assert mutprofile.good_samples(run, 0.05) == {"a"} and mutprofile.good_samples(run, 1) is None and mutprofile.good_samples(run, None) is None
    # N, deletions and uncovered stretches are 'unknown' bases of a genome (minimap2 --cs)
    paf = ["q1\t200\t0\t100\t+\tref\t1000\t0\t100\t90\t100\t60\tcs:Z::10*gn:5-ac:20"]
    gaps = lineages.gap_intervals_from_paf(paf, 150)["q1"]
    assert (11, 11) in gaps and (17, 18) in gaps and (101, 150) in gaps and lineages.is_unknown(gaps, 12) is False and lineages.is_unknown(gaps, 11) and lineages.is_unknown(gaps, 120)
    # Simpson's paradox: later genomes belong to the lineage that carries more inherited mutations; inside each lineage nothing changes over time
    rows = [{"t": 2024.0 + i * 0.05, "apobec_var": 10, "lineage": "X"} for i in range(10)] + [{"t": 2025.0 + i * 0.05, "apobec_var": 30, "lineage": "Y"} for i in range(10)]
    assert outbreak.trend([r["t"] for r in rows], [r["apobec_var"] for r in rows])["slope_per_year"] > 10                     # overall: strong rise
    w = outbreak.within_lineage_trend(rows)
    assert abs(w["slope_per_year"]) < 1e-6 and w["lineages"] == 2                                                             # within lineages: none
    rows2 = [{"t": 2024.0 + i * 0.1, "apobec_var": 10 + 6 * i * 0.1, "lineage": "X"} for i in range(10)]
    assert abs(outbreak.within_lineage_trend(rows2)["slope_per_year"] - 6) < 0.01


def test_origin_classes_opportunities_and_enrichment():
    from genomecurator import insight
    seq = "CCTCAGAGTT"                                   # 1-based: 1 C, 2 C, 3 T, 4 C (after T), 5 A, 6 G (before A), 7 A, 8 G, 9 T, 10 T
    assert insight._classes("C", "T", 4, seq) == "apobec_TC" and insight._classes("G", "A", 6, seq) == "apobec_GA"
    assert insight._classes("C", "T", 2, seq) == "other_CT_GA" and insight._classes("A", "G", 5, seq) == "other_transition" and insight._classes("C", "A", 1, seq) == "transversion"
    opp = insight._opportunities(seq, 1, 10)
    assert opp["C_apobec"] == 1 and opp["C_other"] == 2 and opp["G_apobec"] == 1 and opp["G_other"] == 1 and opp["all"] == 10
    e = insight._enrichment(80, 20, 100, 400)           # 80 % of the C>T/G>A changes in APOBEC contexts that hold 20 % of the C/G sites
    assert e["expected_share"] == 0.2 and e["enrichment"] == 4.0 and e["rate_ratio"] == 16.0 and e["attributable"] == 75       # 80 observed - 20 * 0.25 / ... : background rate 20/400 * 100 sites = 5
    assert insight._enrichment(0, 0, 100, 400) == {} and insight._enrichment(5, 5, 0, 10) == {}


def test_gene_hotspots_use_independent_events_not_inherited_sums():
    import math
    from genomecurator import outbreak
    assert abs(outbreak._pois_sf(0, 3) - 1) < 1e-12 and abs(outbreak._pois_sf(1, 3) - (1 - math.exp(-3))) < 1e-12
    assert abs(outbreak._pois_sf(10, 1.7) - 1.19728e-5) < 1e-9 and outbreak._pois_sf(10, 1.7) < 1e-4 and outbreak._pois_sf(3, 17) > 0.99      # 10 events where 1.7 are expected: a hotspot; 3 where 17 are expected: not


def test_world_map_data_and_country_names(tmp_path, monkeypatch):
    import json
    from genomecurator import geo
    w = json.loads(geo.outlines())
    assert len(w["countries"]) > 150 and "public domain" in w["source"] and all({"iso", "name", "c", "poly"} <= set(c) for c in w["countries"])
    assert geo.resolve("USA")["iso"] == "USA" and geo.resolve("Viet Nam")["iso"] == "VNM" and geo.resolve("Zaire")["iso"] == "COD" and geo.resolve("Cote d'Ivoire")["iso"] == "CIV"
    assert geo.resolve("Democratic Republic of the Congo")["iso"] == "COD" and geo.resolve("Republic of the Congo")["iso"] == "COG"            # the two Congos stay apart
    assert geo.resolve("Singapore")["iso"] == "SGP" and geo.resolve("Hong Kong")["lat"] > 22                                                   # small places the outlines lack
    assert geo.resolve("missing") is None and geo.resolve("") is None and geo.resolve("Atlantis") is None
    _qc_env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    r = TestClient(web.create_app(scheduler=False)).get("/api/world")
    assert r.status_code == 200 and len(r.json()["countries"]) == len(w["countries"])


def test_novel_variant_discovery_statuses_and_known_mutations(monkeypatch):
    from genomecurator import lineages, mutprofile, novel, outbreak
    tip = lambda n: {"name": n}
    ref = "A" * 9 + "G" + "A" * 9 + "T" + "A" * 9 + "G" + "A" * 9000
    tree = {"name": "root", "children": [
        {"name": "X", "branch_attrs": {"mutations": {"nuc": ["G10A", "T20C"]}}, "children": [tip("x1"), tip("x2"), tip("x3"), tip("x4")]},        # G10A is a defining mutation: only T20C is new
        {"name": "Y", "branch_attrs": {"mutations": {"nuc": ["G30A"]}}, "children": [tip("y1"), tip("y2"), tip("y3")]},
        {"name": "Z", "branch_attrs": {"mutations": {"nuc": ["G10A"]}}, "children": [tip("z1"), tip("z2"), tip("z3")]}]}                       # only a known mutation: no candidate
    monkeypatch.setattr(outbreak, "_tree", lambda run: (tree, ref, {"name": "t", "build": "x-no-such-build"}))
    monkeypatch.setattr(outbreak, "_tree_reference", lambda ns: ("R", ref, []))
    monkeypatch.setattr(novel, "_known_in_tree_coordinates", lambda fa: ({(10, "A")}, True))
    dt = {"x1": "2025-01-05", "x2": "2025-02-05", "x3": "2025-03-05", "x4": "2025-04-05", "y1": "2026-08-01", "y2": "2026-08-20", "y3": "2026-09-10", "z1": "2026-09-01", "z2": "2026-09-02", "z3": "2026-09-03"}
    gs = [{"id": i, "lineage": "A", "outbreak": "o", "clade": "Ib", "country": "U" if i[0] != "y" else "K", "date": d} for i, d in dt.items()]
    monkeypatch.setattr(lineages, "_genomes", lambda run: gs)
    monkeypatch.setattr(mutprofile, "good_samples", lambda run, mn: None)
    d = novel.discover({"dataset_id": 1}, window_days=120, min_size=3)
    byfirst = {c["mutations"][0]["label"]: c for c in d["candidates"]}
    assert set(byfirst) == {"T20C", "G30A"}                                           # the defining G10A (on X and Z) makes neither a candidate
    assert byfirst["G30A"]["status"] == "new" and byfirst["G30A"]["countries"] == {"K": 3} and byfirst["G30A"]["n_recent"] == 3
    assert byfirst["T20C"]["status"] == "steady" and byfirst["T20C"]["n"] == 4 and d["counts"] == {"new": 1, "steady": 1}
    assert d["candidates"][0]["status"] == "new" and d["n_novel_mutations"] == 2 and d["recent_novel_mutations"] == 1                     # new ones are listed first
    assert d["curve"][-1]["total"] == 2 and {c["month"]: c["new"] for c in d["curve"]}["2025-01"] == 1
    assert novel.discover({"dataset_id": 1}, min_size=5)["candidates"] == []          # a larger minimum size drops them


def test_mutation_groups_and_pca_of_genomes(monkeypatch):
    import numpy as np
    from genomecurator import cooccur, insight, lineages, mutprofile, novel
    d = np.array([[0, .02, .9, .9], [.02, 0, .9, .9], [.9, .9, 0, .1], [.9, .9, .1, 0]])
    assert sorted(sorted(c) for c in cooccur._upgma(d, 0.15)) == [[0, 1], [2, 3]] and len(cooccur._upgma(d, 0.0)) == 4 and len(cooccur._upgma(d, 0.95)) == 1
    hdr = ["sample", "pos", "ref", "alt", "gene", "aa_change", "effect", "profile_of_interest"]
    rows = []
    def add(pos, samples, apo="other"):
        for s in samples:
            rows.append({"sample": s, "pos": str(pos), "ref": "G", "alt": "A", "gene": "G1", "aa_change": "K2R", "effect": "nonsynonymous", "profile_of_interest": apo})
    A = [f"a{i}" for i in range(6)]; B = [f"b{i}" for i in range(6)]
    add(100, A); add(200, A); add(300, A[:5] + ["b0"])             # 100 and 200 have identical carriers, 300 almost (one genome swapped)
    add(500, B, "GA>AA"); add(600, B)                             # a second group
    add(900, ["a0", "b0", "a1"])                                  # a mutation that does not follow the groups
    monkeypatch.setattr(insight, "_mutations", lambda run, mn=None: rows)
    monkeypatch.setattr(lineages, "_genomes", lambda run: [{"id": s, "lineage": "X" if s[0] == "a" else "Y", "outbreak": "o", "clade": "", "country": "U", "date": "2025-01-01"} for s in A + B])
    monkeypatch.setattr(mutprofile, "summary", lambda run: {"reference": "R"})
    monkeypatch.setattr(novel, "discover", lambda *a, **k: {"available": False})
    r = cooccur.analyse({"dataset_id": 1}, max_n=None, min_carriers=3, merge=0.3)
    assert r["available"] and r["n_genomes"] == 12 and r["n_sites"] == 6
    groups = {frozenset(m["label"] for m in g["mutations"]) for g in r["groups"]}
    assert frozenset({"G100A", "G200A", "G300A"}) in groups and frozenset({"G500A", "G600A"}) in groups and frozenset({"G900A"}) in groups
    assert r["explained"][0] > 0.4 and r["lineage_eta2"][0] > 0.8                                    # the first axis is the lineage split
    x = {g["id"]: g["pc"][0] for g in r["genomes"]}
    assert all(x[a] * x["b5"] < 0 for a in A[1:5])                                                    # lineages X and Y lie on opposite sides of PC1
    strict = cooccur.analyse({"dataset_id": 1}, max_n=None, min_carriers=3, merge=0.0)
    assert frozenset({"G100A", "G200A"}) in {frozenset(m["label"] for m in g["mutations"]) for g in strict["groups"]} and len(strict["groups"]) > len(r["groups"])
    assert cooccur.analyse({"dataset_id": 1}, max_n=None, min_carriers=50)["available"] is False


def test_lineage_heatmap_frequencies_defining_marks_and_gvf(tmp_path, monkeypatch):
    import io, zipfile
    from genomecurator import heatmap, insight, lineages, mutprofile
    ref = tmp_path / "r.fasta"
    ref.write_text(">R\n" + "ACGT" * 60 + "\n", encoding="utf-8")
    monkeypatch.setattr(mutprofile, "reference_by_name", lambda acc: (ref, None))
    monkeypatch.setattr(mutprofile, "summary", lambda run: {"reference": "R"})
    monkeypatch.setattr(mutprofile, "good_samples", lambda run, mn: None)
    monkeypatch.setattr(mutprofile, "run_itr", lambda run: None)
    row = lambda s, p, a="A", apo="other": {"sample": s, "pos": str(p), "ref": "G", "alt": a, "gene": "G1", "aa_change": "K2R", "effect": "nonsynonymous", "profile_of_interest": apo}
    X = [f"x{i}" for i in range(4)]; Y = [f"y{i}" for i in range(4)]
    rows = [row(s, 10) for s in X] + [row(s, 20, apo="GA>AA") for s in X[:3]] + [row(s, 10) for s in Y[:1]] + [row(s, 30) for s in Y]
    monkeypatch.setattr(insight, "_mutations", lambda run, mn=None: rows)
    gs = [{"id": s, "lineage": "X" if s[0] == "x" else "Y", "outbreak": "o", "clade": "Ib", "country": "U", "date": "2025-01-01"} for s in X + Y]
    monkeypatch.setattr(lineages, "_genomes", lambda run: gs)
    monkeypatch.setattr(lineages, "status", lambda: {"ok": True, "commit": "t"})
    monkeypatch.setattr(lineages, "defining_in_run_coordinates", lambda run, fa: None)
    monkeypatch.setattr(lineages, "defining_snps", lambda ob, lin: [(10, "A", "X")] if lin == "X" else [])
    h = heatmap.build({"dataset_id": 1}, "lineage", min_n=2, max_n=None, min_freq=0.0, min_carriers=1)
    assert {r["id"]: r["n"] for r in h["rows"]} == {"X": 4, "Y": 4} and [c["label"] for c in h["columns"]] == ["G10A", "G20A", "G30A"]
    cell = {(h["rows"][r]["id"], h["columns"][c]["label"]): (f, n, d) for r, c, f, n, d in h["cells"]}
    assert cell[("X", "G10A")] == (100.0, 4, 1) and cell[("Y", "G10A")] == (25.0, 1, 0)                  # frequency within the lineage; defining only for the lineage it defines
    assert cell[("X", "G20A")][0] == 75.0 and ("Y", "G20A") not in cell and cell[("Y", "G30A")][0] == 100.0
    assert h["columns"][1]["apobec"] and h["columns"][0]["defined_by"] == "X" and h["rows"][[r["id"] for r in h["rows"]].index("X")]["n_defining"] == 1
    assert heatmap.build({"dataset_id": 1}, "lineage", min_n=5, max_n=None)["rows"] == []              # groups under the minimum size are left out
    text = heatmap.gvf({"dataset_id": 1}, "X", "lineage", None)
    assert text.startswith("##gff-version 3") and "##gvf-version 1.10" in text and "sample_group=X" in text
    body = [l for l in text.splitlines() if not l.startswith("#")]
    assert len(body) == 2 and "ao=4;dp=4;sample_size=4" in body[0] and "clade_defining=True" in body[0] and "apobec3_type=True" in body[1]
    z = zipfile.ZipFile(io.BytesIO(heatmap.gvf({"dataset_id": 1}, "", "lineage", None)))
    assert sorted(z.namelist()) == ["X_annotated.gvf", "Y_annotated.gvf"]
    with __import__("pytest").raises(ValueError):
        heatmap.gvf({"dataset_id": 1}, "nope", "lineage", None)


def test_menu_groups_hub_and_result_tab_strip(tmp_path, monkeypatch):
    _qc_env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    c = TestClient(web.create_app(scheduler=False))
    home = c.get("/").text
    assert home.count('<nav class="subnav">') == 1 and "변경 검토" in home and "정합성" in home                 # the active group's sub-tabs are shown
    for url, group, tab in (("/genomes", "유전체", "유전체 표"), ("/datasets", "데이터셋", "내보내기"), ("/hub", "분석", "분석 허브"), ("/nextstrain", "분석", "실행·이력"),
                            ("/pathoplexus", "제출·공유", "Superset"), ("/jobs", "시스템", "작업"), ("/analysis", "제출·공유", "Pathoplexus 제출")):
        t = c.get(url).text
        sub = t[t.index('<nav class="subnav">'):t.index("</nav>")]
        assert f'<span class="sg">{group}</span>' in sub and tab in sub, url
    assert c.get("/hub").status_code == 200 and 'id="jbadge"' in c.get("/hub").text                          # the jobs badge survives on every page
    assert c.get("/analysis/999/snipit").status_code == 404                                                  # unknown runs are still refused
    # the result pages carry the shared tab strip with every view
    import pathlib
    for name in ("snipit", "heatmap", "pca", "spread", "insight", "outbreak", "novel"):
        src = (pathlib.Path(web.__file__).parent / "templates" / f"{name}.html").read_text(encoding="utf-8")
        assert "{% set view='" + name + "' %}{% include \"_anav.html\" %}" in src, name
    anav = (pathlib.Path(web.__file__).parent / "templates" / "_anav.html").read_text(encoding="utf-8")
    assert all(f"('{v}'," in anav for v in ("snipit", "heatmap", "pca", "spread", "insight", "outbreak", "novel"))


def test_page_scripts_have_no_syntax_errors_and_help_is_collapsible(tmp_path, monkeypatch):
    """Every inline script of the analysis pages must parse (a raw newline inside a JS string once broke a page silently); needs node, skipped without it."""
    import re, shutil, subprocess
    node = shutil.which("node")
    if not node:
        import pytest
        pytest.skip("node is not installed")
    from pathlib import Path
    from genomecurator import web
    tdir = Path(web.__file__).parent / "templates"
    for name in ("snipit", "heatmap", "pca", "spread", "insight", "outbreak", "novel", "hub", "vocab", "qc", "concordance", "base"):
        html = (tdir / f"{name}.html").read_text(encoding="utf-8")
        for k, body in enumerate(re.findall(r"<script>(.*?)</script>", html, re.S)):
            if "{%" in body or "{{" in body:                    # Jinja inside: render-free check of the plain parts only
                body = re.sub(r"\{\{.*?\}\}", "0", body, flags=re.S)
                body = re.sub(r"\{%.*?%\}", "", body, flags=re.S)
            f = tmp_path / f"{name}_{k}.js"
            f.write_text(body, encoding="utf-8")
            r = subprocess.run([node, "--check", str(f)], capture_output=True, text=True)
            assert r.returncode == 0, f"{name} script {k}: {r.stderr[:300]}"
    base = (tdir / "base.html").read_text(encoding="utf-8")
    assert "function gcSyncN" in base.split("</head>")[0] and "function gcHelp" in base                           # the helpers the page scripts call are defined before them


def test_snvs_from_a_squirrel_alignment_and_the_aligner_choice(tmp_path, monkeypatch):
    import gzip
    from genomecurator import aligned, lineages, mutprofile, squirrel
    # target reference T (14 bp) and Squirrel's own reference A (the same sequence with 2 extra bases at the start, so columns are shifted by 2)
    T = "ACGTACGTACGTAC"
    A = "GG" + T
    tfa, afa = tmp_path / "t.fasta", tmp_path / "a.fasta"
    tfa.write_text(">T\n" + T + "\n", encoding="utf-8")
    afa.write_text(">A\n" + A + "\n", encoding="utf-8")
    monkeypatch.setattr(aligned, "_liftover", lambda q, t: {i: i - 2 for i in range(3, len(A) + 1)})          # column i of A is position i-2 of T
    rows = {"s1": A, "s2": A[:6] + "G" + A[7:], "s3": A[:8] + "NN" + A[10:], "s4": A[:3] + "T" + A[4:]}      # s2: T pos 5 (A) -> G; s3: two N's; s4: T pos 2 (C) -> T
    aln = tmp_path / "aln.fasta.gz"
    with gzip.open(aln, "wt") as f:
        for k, v in rows.items():
            f.write(f">{k} description\n{v}\n")
    r = aligned.snvs_vs_target(aln, afa, tfa, ids={"s1", "s2", "s3", "s4"})
    got = {(n, p, ref, alt) for n, p, ref, alt in r["snvs"]}
    assert ("s2", 5, "A", "G") in got and ("s4", 2, "C", "T") in got and not [x for x in got if x[0] in ("s1", "s3")]       # identical genome: none; N is never a difference
    assert any(a <= 7 <= b for a, b in r["gaps"]["s3"]) and any(a <= 8 <= b for a, b in r["gaps"]["s3"])                    # its N columns are unknown ...
    assert r["gaps"]["s1"] == []                                                                                         # a clean genome with every reference position covered has none
    assert r["n"] == 4 and sorted(aligned.read_alignment(aln))[0][0] == "s1"
    assert [x for x in aligned.snvs_vs_target(aln, afa, tfa, ids={"s2"})["snvs"]] == [("s2", 5, "A", "G")]                  # only the requested genomes
    with __import__("pytest").raises(ValueError):                                                                           # a row of the wrong length is refused
        bad = tmp_path / "bad.fasta"
        bad.write_text(">x\nACGT\n", encoding="utf-8")
        aligned.snvs_vs_target(bad, afa, tfa)
    # which alignment is used
    monkeypatch.setattr(squirrel, "status", lambda: {"ok": True})
    monkeypatch.setattr(squirrel, "reference_fasta", lambda acc: afa)
    from genomecurator import pathoplexus
    monkeypatch.setattr(pathoplexus, "organism_for", lambda w: "mpox" if w == 1 else "")
    ds = {"watch_id": 1, "build": "clade-i", "squirrel_json": ""}
    assert mutprofile.aligner_choice(ds) == "squirrel" and squirrel.alignment_reference(ds) == "PX667572"
    assert mutprofile.aligner_choice({**ds, "watch_id": 2}) == "minimap2"                                                    # other viruses keep minimap2
    assert mutprofile.aligner_choice({**ds, "mutprofile_json": '{"aligner": "minimap2"}'}) == "minimap2"
    assert mutprofile.settings_sha({**ds, "mutprofile_json": '{"aligner": "minimap2"}'}) != mutprofile.settings_sha({**ds, "mutprofile_json": '{"aligner": "auto"}'})   # another aligner is another analysis
    assert mutprofile.aligner_choice({"watch_id": 2, "build": "", "mutprofile_json": '{"aligner": "squirrel"}', "squirrel_json": '{"clade": "cladeii"}'}) == "squirrel"
    assert squirrel.alignment_reference({"build": "", "squirrel_json": '{"clade": "cladeiib"}'}) == "NC_063383" and squirrel.alignment_reference({"build": "", "squirrel_json": '{"clade": "split"}'}) is None
    monkeypatch.setattr(squirrel, "reference_fasta", lambda acc: None)
    with __import__("pytest").raises(RuntimeError):
        mutprofile.aligner_choice({**ds, "mutprofile_json": '{"aligner": "squirrel"}'})
    assert mutprofile.parse_settings({"aligner": "squirrel"}) == {"aligner": "squirrel"}
    with __import__("pytest").raises(ValueError):
        mutprofile.parse_settings({"aligner": "mafft"})


def test_vocab_standardises_values_without_guessing(tmp_path, monkeypatch):
    from genomecurator import config, vocab
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    vocab._CACHE.clear()
    S = vocab.standardize
    assert vocab.key("Côte d'Ivoire") == vocab.key("Cote d’Ivoire") == "cote d ivoire" and vocab.key("H. sapiens") == "h sapiens"
    # host: dictionary synonyms, the typo we saw, non-hosts kept and flagged, unknown kept and marked
    h = S("host_name", "Homosapiens")
    assert (h.value, h.status, h.rule) == ("Homo sapiens", "mapped", "host.synonym")
    assert S("host_name", "Homo sapiens").status == "same" and S("host_name", " homo   sapiens ").value == "Homo sapiens"
    assert S("host_name", "Patent").value == "Homo sapiens" and S("host_name", "cynomolgus monkey").value == "Macaca fascicularis"
    v = S("host_name", "Vero e6")
    assert v.status == "flagged" and v.value == "Vero E6" and "cell line" in v.rule
    u = S("host_name", "Mus imaginarius")
    assert u.status == "unmapped" and u.value == "Mus imaginarius"                                       # never replaced, never dropped
    assert S("host_name", "missing").status == "empty" and S("host_name", "NA").value == ""
    # country
    assert S("isolation_country", "Zaire").value == "Democratic Republic of the Congo" and S("isolation_country", "Viet Nam").status == "same" and S("isolation_country", "Vietnam").value == "Viet Nam"
    assert S("isolation_country", "United States").value == "USA" and S("isolation_country", "USA: NY, New York").value == "USA" and S("isolation_country", "Atlantis").status == "unmapped"
    assert S("isolation_country", "Congo").value == "Republic of the Congo" and S("isolation_country", "DRC").value == "Democratic Republic of the Congo"        # the two Congos stay apart
    # clade: only for mpox
    assert S("clade", "Clade Ib", "mpox").value == "Ib" and S("clade", "Ib", "mpox").status == "same" and S("clade", "clade 2b", "mpox").value == "IIb"
    assert S("clade", "Ib/IIb", "mpox").status == "flagged" and S("clade", "II", "mpox").status == "flagged" and S("clade", "xyz", "mpox").status == "unmapped"
    assert S("clade", "whatever", "*").status == "nodict" and S("clade", "whatever", "*").value == "whatever"        # no clade dictionary for other viruses: untouched
    # dates
    D = vocab.standardize_date
    assert D("2025-05-01").status == "same" and D("2025-05").status == "same" and D("2025").status == "same"
    assert (D("2025-05-01/2025-05-28").value, D("2025-05-01/2025-05-28").status) == ("2025-05", "normalized") and D("2024-01-01/2024-12-09").value == "2024"
    assert D("2023-12-01/2024-02-01").status == "flagged" and D("2023-12-01/2024-02-01").value == "2023-12-01/2024-02-01"          # across years: not reduced to one date
    assert D("15-May-2025").value == "2025-05-15" and D("May 2025").value == "2025-05" and D("2025/05/15").value == "2025-05-15" and D("20250515").value == "2025-05-15" and D("2025-05-XX").value == "2025-05"
    assert D("13/05/2025").value == "2025-05-13" and D("05/13/2025").value == "2025-05-13"
    assert D("05/06/2025").status == "flagged" and D("05/06/2025").value == "05/06/2025"                                  # day / month order unknown: left alone
    assert D("2025-02-31").status == "flagged" and D("null").status == "empty" and D("soon").status == "flagged"
    # a user's file overrides and extends the shipped dictionary
    (tmp_path / "vocab").mkdir()
    (tmp_path / "vocab" / "host.tsv").write_text("term\tsynonyms\ttaxon_id\tflag\tnote\nMus imaginarius\tfake mouse\t1\t\t\nHomo sapiens\thuman being\t9606\t\t\n", encoding="utf-8")
    vocab._CACHE.clear()
    assert S("host_name", "fake mouse").value == "Mus imaginarius" and S("host_name", "Mus imaginarius").status == "same" and S("host_name", "human being").value == "Homo sapiens" and S("host_name", "Homosapiens").value == "Homo sapiens"
    assert vocab.scope_for("Orthopoxvirus monkeypox", "") == "mpox" and vocab.scope_for("Henipavirus nipahense", "") == "*"


def test_check_standardises_before_sources_are_compared_and_keeps_the_raw_value(tmp_path, monkeypatch):
    import json
    from genomecurator import config, curate, sources, vocab
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "g.db")
    monkeypatch.setattr(config, "MIN_COVERAGE_PCT", 0.0)
    vocab._CACHE.clear()
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon,taxon_id) VALUES(1,'Orthopoxvirus monkeypox','Orthopoxvirus monkeypox',10244)")
    con.commit()
    con.close()
    rec = lambda **k: {"accession_version": 1, "genome_length": 1000, "genome_status": "Complete", **k}
    ncbi = {"A1": rec(host_name="Homosapiens", isolation_country="Zaire", collection_date="2025-05-01/2025-05-28", clade="Clade Ib"), "A2": rec(host_name="Vero e6", isolation_country="Atlantis")}
    ppx = {"A1": {"host_name": "Homo sapiens", "isolation_country": "Democratic Republic of the Congo", "clade": "Ib"}}
    monkeypatch.setattr(sources, "ncbi_species_scope", lambda tid: None)
    monkeypatch.setattr(sources, "pathoplexus_resolve", lambda *a, **k: None)
    monkeypatch.setattr(sources, "fetch_all", lambda *a, **k: ({"ncbi": json.loads(json.dumps(ncbi)), "bvbrc": {}, "pathoplexus": json.loads(json.dumps(ppx))}, []))
    cs = curate.check("Orthopoxvirus monkeypox")
    con = connect()
    new = {r["accession"]: json.loads(r["record_json"] if "record_json" in r.keys() else "{}") for r in con.execute("SELECT * FROM changeset_item WHERE changeset_id=?", (cs,))}
    std = {(r["accession"], r["source"], r["field"]): dict(r) for r in con.execute("SELECT * FROM standardization WHERE watch_id=1")}
    con.close()
    r = std[("A1", "ncbi", "host_name")]
    assert (r["raw"], r["standard"], r["status"]) == ("Homosapiens", "Homo sapiens", "mapped")                          # the raw value is kept, with the reason
    assert std[("A1", "ncbi", "isolation_country")]["raw"] == "Zaire" and std[("A1", "ncbi", "collection_date")]["standard"] == "2025-05"
    assert std[("A1", "ncbi", "clade")]["standard"] == "Ib" and ("A1", "pathoplexus", "host_name") not in std       # already standard values are not recorded
    assert std[("A2", "ncbi", "host_name")]["status"] == "flagged" and std[("A2", "ncbi", "isolation_country")]["status"] == "unmapped"
    # the sources now agree: NCBI says Zaire / Homosapiens, Pathoplexus the standard names - no conflict and no notation difference is reported for A1
    con = connect()
    rows = [dict(x) for x in con.execute("SELECT * FROM source_conflict WHERE watch_id=1 AND accession='A1' AND field IN ('host_name','isolation_country','clade')")]
    con.close()
    assert rows == []
    # the curated record carries the standard values
    con = connect()
    item = con.execute("SELECT * FROM changeset_item WHERE changeset_id=? AND accession='A1'", (cs,)).fetchone()
    blob = " ".join(str(item[k]) for k in item.keys() if k.endswith("_json"))
    con.close()
    assert "Homo sapiens" in blob and "Democratic Republic of the Congo" in blob and "Homosapiens" not in blob and "Zaire" not in blob


def test_squirrel_alignment_reprojected_on_the_workflow_reference(tmp_path, monkeypatch):
    import gzip
    from genomecurator import aligned, lineages
    T = "ACGTACGTAC"                                   # the workflow's reference (10 bp)
    A = "GG" + T[:5] + "TT" + T[5:]                    # Squirrel's reference: 2 extra bases in front and 2 in the middle that the workflow reference does not have
    tfa, afa = tmp_path / "t.fasta", tmp_path / "a.fasta"
    tfa.write_text(">T\n" + T + "\n", encoding="utf-8")
    afa.write_text(">A\n" + A + "\n", encoding="utf-8")
    cmap = {**{i: i - 2 for i in range(3, 8)}, **{i: i - 4 for i in range(10, 15)}}                  # column of A -> position of T (A cols 8-9 have no counterpart)
    monkeypatch.setattr(aligned, "_liftover", lambda q, t: cmap)
    aln = tmp_path / "aln.fasta.gz"
    seq1 = A[:5] + "C" + A[6:]                                                                        # a change at A col 6 (T pos 4)
    with gzip.open(aln, "wt") as f:
        f.write(f">g1 x\n{A}\n>g2\n{seq1}\n>g3\n{A[:10]}NNNN\n")
    got = dict(aligned.reproject(aln, afa, tfa))
    assert set(got) == {"g1", "g2", "g3"} and all(len(v) == len(T) for v in got.values())             # exactly the reference's columns, like Nextclade's --output-fasta
    assert got["g1"] == T and got["g2"][3] == "C" and got["g2"][:3] + got["g2"][4:] == T[:3] + T[4:]
    assert got["g3"] == T[:6] + "NNNN" and set(dict(aligned.reproject(aln, afa, tfa, ids={"g2"}))) == {"g2"}
    cmap.clear()
    cmap.update({i: i for i in range(1, 4)})                                                          # positions the column map does not reach are N for every genome
    assert dict(aligned.reproject(aln, afa, tfa))["g1"] == A[:3] + "N" * 7


def test_vocab_editing_decisions_land_in_the_user_dictionary(tmp_path, monkeypatch):
    import pytest
    from fastapi.testclient import TestClient
    from genomecurator import config, vocab, web
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    vocab._CACHE.clear()
    S = vocab.standardize
    assert S("host_name", "Mus fantasticus").status == "unmapped"
    vocab.add_term("host_name", "Mus fantasticus", ["fantasy mouse"], extra={"taxon_id": "1"})
    assert S("host_name", "fantasy  Mouse").value == "Mus fantasticus"
    # a writing cannot belong to two terms
    with pytest.raises(ValueError):
        vocab.add_term("host_name", "Other host", ["fantasy mouse"])
    with pytest.raises(ValueError):
        vocab.add_synonym("host_name", "x", "No such term")
    vocab.add_synonym("host_name", "Humans", "Homo sapiens")
    assert S("host_name", "humans").value == "Homo sapiens"
    # a flagged shipped term accepted as standard: its own spellings lose the flag, and removing the user's row restores the shipped state
    assert S("host_name", "Vero E6").status == "flagged"
    vocab.add_term("host_name", "Vero E6", [], "", "ok")
    assert S("host_name", "Vero E6").status != "flagged"
    assert vocab.remove_entry("host_name", "Vero E6") and S("host_name", "Vero E6").status == "flagged" and not vocab.remove_entry("host_name", "Vero E6")
    # dates: only a valid date is accepted, and the decision is applied
    assert vocab.standardize_date("05/06/2024").status in ("flagged", "unmapped")
    with pytest.raises(ValueError):
        vocab.set_date("05/06/2024", "2024-13-40")
    vocab.set_date("05/06/2024", "2024-06-05")
    assert vocab.standardize_date("05/06/2024").value == "2024-06-05"
    assert [m["term"] for m in vocab.suggest("host_name", "Homo sapien")][:1] == ["Homo sapiens"] and vocab.suggest("host_name", "zzzzzz") == []
    assert any(e["term"] == "Mus fantasticus" and e["user"] for e in vocab.search_terms("host_name", "fantasy"))
    # the same through the web API
    c = TestClient(web.create_app(scheduler=False))
    r = c.post("/api/vocab/term", json={"field": "host_name", "term": "Rattus fantasticus", "synonyms": ["fantasy rat"]})
    assert r.status_code == 200 and S("host_name", "fantasy rat").value == "Rattus fantasticus"
    assert c.post("/api/vocab/term", json={"field": "host_name", "term": "Dup", "synonyms": ["fantasy rat"]}).status_code == 400
    assert c.post("/api/vocab/date", json={"raw": "7/8/2023", "standard": "2023-08-07"}).status_code == 200
    assert any(e["term"] == "Rattus fantasticus" for e in c.get("/api/vocab/terms", params={"field": "host_name", "q": "rattus"}).json())
    assert c.post("/api/vocab/remove", json={"field": "host_name", "term": "Rattus fantasticus"}).json() == {"removed": True}
    assert c.get("/api/vocab/terms", params={"field": "nonsense"}).status_code == 400


def test_consensus_mask_marks_only_columns_where_aligners_disagree(tmp_path):
    from genomecurator import aligned
    ref = "ACGTACGTACGTACGTACGT"
    (tmp_path / "ref.fa").write_text(f">r\n{ref}\n")
    def mut(s, i, b):
        return s[:i] + b + s[i + 1:]
    sq = {"g1": mut(ref, 4, "T"), "g2": mut(ref, 4, "T"), "g3": mut(ref, 4, "T"), "g4": mut(mut(ref, 9, "N"), 12, "C")}
    nc = {"g1": ref, "g2": ref, "g3": ref, "g4": mut(ref, 9, "C")}                                  # column 5: three genomes differ; 10: N vs base = unknown, not a disagreement; 13: only g4 differs
    (tmp_path / "sq.fa").write_text("".join(f">{k}\n{v}\n" for k, v in sq.items()))
    (tmp_path / "nc.fa").write_text("".join(f">{k}\n{v}\n" for k, v in nc.items()))
    d = aligned.disagreement(tmp_path / "sq.fa", tmp_path / "ref.fa", tmp_path / "nc.fa", tmp_path / "ref.fa")
    assert d["columns"] == {5: 3, 13: 1} and d["n"] == 4
    assert aligned.consensus_mask(d, 3) == [(5, 5)] and aligned.consensus_mask(d, 1) == [(5, 5), (13, 13)]
    assert aligned.mask_bed([(5, 5), (13, 14)], "r") == "r\t4\t5\nr\t12\t14\n"


def test_variantcall_vcf_parsing_and_concordance(tmp_path):
    from genomecurator import variantcall
    v = tmp_path / "x.vcf"
    v.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nr\t10\t.\tA\tG\t50\tPASS\tDP=17;AF=1\nr\t20\t.\tAT\tA\t50\tPASS\tDP=9\nr\t30\t.\tC\tT\t50\tPASS\tDP=20;AF=0.4\n")
    p = variantcall.parse_vcf(v)
    assert set(p) == {(10, "A", "G"), (30, "C", "T")} and p[(30, "C", "T")] == {"dp": 20, "af": 0.4}
    c = variantcall.concordance({"a": {"g": {(10, "A", "G"), (30, "C", "T")}}, "b": {"g": {(10, "A", "G")}}, "c": {"h": set()}})
    by = {r["pos"]: r for r in c["rows"]}
    assert by[10]["n"] == 2 and by[10]["of"] == 2 and by[30]["by"] == ["a"] and by[30]["of"] == 2        # genome g was judged by a and b only
    assert c["summary"]["all_agree"] == 1 and c["summary"]["only_one"] == 1


def test_concordance_for_run_lists_only_sites_where_callers_differ(tmp_path, monkeypatch):
    from genomecurator import mutprofile, variantcall
    rd = tmp_path / "mutprofile"
    rd.mkdir()
    (rd / "mutations.tsv").write_text("sample\tpos\tref\talt\ng1\t10\tA\tG\ng1\t20\tC\tT\ng2\t10\tA\tG\n")
    (rd / "variantcall.tsv").write_text("genome\tpos\tref\talt\tsnippy\tlofreq\tdp\taf\ng1\t10\tA\tG\t1\t1\t17\t1\ng1\t20\tC\tT\t0\t1\t17\t1\ng2\t10\tA\tG\t1\t1\t17\t1\ng2\t30\tG\tA\t1\t0\t\t\n")
    (rd / "variantcall.genomes.txt").write_text("g1\ng2\n")
    monkeypatch.setattr(mutprofile, "run_itr", lambda run: None)
    monkeypatch.setattr(mutprofile, "good_samples", lambda run, mx: None)
    c = variantcall.concordance_for_run({"run_dir": str(tmp_path)}, None)
    assert c["summary"]["variants"] == 4 and c["summary"]["all_agree"] == 2
    assert [(s["pos"], s["mutprofile"], s["snippy"], s["lofreq"]) for s in c["sites"]] == [(20, 1, 0, 1), (30, 0, 1, 0)]
    assert not variantcall.available({"run_dir": str(tmp_path / "x")})


def test_vocab_platform_and_specimen_fold_into_standard_terms(tmp_path, monkeypatch):
    from genomecurator import config, vocab
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    vocab._CACHE.clear()
    S = vocab.standardize
    # the model stays in sequencing_platform (spelling made standard); the technology goes to the derived column `platform`
    assert S("sequencing_platform", "Illumina NovaSeq 6000").status == "same" and S("sequencing_platform", "novaseq").value == "Illumina NovaSeq 6000"
    assert S("sequencing_platform", "Oxford Nanopore Technologies GridION").value == "Oxford Nanopore GridION" and S("sequencing_platform", "MinION").value == "Oxford Nanopore MinION"
    assert S("sequencing_platform", "Nanopore").value == "Oxford Nanopore"
    r = S("sequencing_platform", "Illumina MiSeq; Oxford Nanopore Technologies GridION")
    assert (r.value, r.status) == ("Illumina MiSeq; Oxford Nanopore GridION", "mapped")
    assert vocab.platform_of(r.value) == "Illumina; Oxford Nanopore" and vocab.platform_of("Illumina NextSeq 2000") == "Illumina" and vocab.platform_of("Quantum") == ""
    assert S("sequencing_platform", "Nanopore+Illumina").value == "Oxford Nanopore; Illumina"
    assert S("sequencing_platform", "WGS").status == "flagged"
    u = S("sequencing_platform", "Illumina; Quantum Sequencer")                                            # one unknown part: kept as it is, for a person to decide
    assert (u.value, u.status) == ("Illumina; Quantum Sequencer", "unmapped")
    recs = {"ncbi": {"A1": {"sequencing_platform": "MinION"}, "A2": {"sequencing_platform": "Quantum"}}}
    vocab.apply_to_sources(recs)
    assert recs["ncbi"]["A1"] == {"sequencing_platform": "Oxford Nanopore MinION", "platform": "Oxford Nanopore"} and "platform" not in recs["ncbi"]["A2"]
    assert S("isolation_source", "lesion swab").value == "Lesion swab" and S("isolation_source", "Skin lesion:Crust").value == "Skin lesion (crust/scab)"
    assert S("isolation_source", "clinical").status == "unmapped"                                         # too vague to guess a specimen type



def test_squirrel_qc_is_compared_with_nextclade_bad_only(tmp_path, monkeypatch):
    import sqlite3
    from genomecurator import config, db, qc
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    con = db.connect()
    con.execute("INSERT INTO watchlist(id,name,taxon) VALUES(1,'V','t')")
    for i, (acc, ppx, st) in enumerate([("A1", "PP_1.1", "bad"), ("A2", "PP_2.1", "bad"), ("A3", "PP_3.1", "mediocre"), ("A4", "", "good")]):
        con.execute("INSERT INTO genome(genome_id,watch_id,accession,pathoplexus_accession) VALUES(?,?,?,?)", (f"g{i}", 1, acc, ppx))
        con.execute("INSERT INTO qc_result(accession,tool,watch_id,status,score) VALUES(?,?,?,?,0)", (acc, "nextclade", 1, st))
    rd = tmp_path / "run" / "squirrel"
    rd.mkdir(parents=True)
    (rd / "suggested_to_exclude.csv").write_text("name,note\nPP_1,N content is 0.3\nA4,N content is 0.4\n")
    (rd / "alignment.suggested_mask.csv").write_text("Name,Minimum,Maximum,Length,present_in,note\n5,5,5,1,PP_1|PP_3,clustered_snps\n")
    (rd / "alignment.aln.fasta.gz").write_bytes(__import__("gzip").compress(b">PP_1\nACGT\n>PP_2\nACGT\n>PP_3\nACGT\n>A4\nACGT\n"))
    con.execute("INSERT INTO analysis_run(id,tool,status,watch_id,dataset_id,run_dir,finished,name,started) VALUES(1,'squirrel','done',1,1,?,'2026-01-01T00:00:00','n','2026-01-01T00:00:00')", (str(tmp_path / "run"),))
    con.commit(); con.close()
    q = qc.squirrel_qc(1)
    assert q["judged"] == 4 and q["flagged"] == 2 and q["both"] == 1 and q["only_squirrel"] == 1 and q["only_nextclade"] == 1 and q["neither"] == 1
    assert q["nextclade_mediocre"] == 1 and q["site_notes"] == {"clustered_snps": 2} and q["genomes_with_site_notes"] == 2
    assert qc.squirrel_qc(99) == {"available": False}


def test_support_badges_say_what_is_built_for_which_virus(tmp_path, monkeypatch):
    from genomecurator import support
    m = {s["area"]: s["level"] for s in support.for_virus("Orthopoxvirus monkeypox")}
    o = {s["area"]: s["level"] for s in support.for_virus("Henipavirus nipahense")}
    assert m == {"collect": "done", "vocab": "done", "qc": "done", "analysis": "done"}
    assert o["analysis"] == "none" and o["vocab"] == "wip" and o["qc"] == "wip"
    from fastapi.testclient import TestClient
    from genomecurator import config, web
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    c = TestClient(web.create_app(scheduler=False))
    assert c.get("/").status_code == 200 and c.get("/hub").status_code == 200


def test_variantflow_parity_and_status(tmp_path, monkeypatch):
    from genomecurator import variantflow
    d = tmp_path / "mutprofile"
    f = tmp_path / "variantflow"
    d.mkdir(); f.mkdir()
    head = "genome\tpos\tref\talt\tsnippy\tlofreq\tdp\taf\n"
    (d / "variantcall.tsv").write_text(head + "g1\t10\tA\tG\t1\t1\t17\t1.0\ng1\t20\tC\tT\t0\t1\t17\t1.0\n")
    run = {"run_dir": str(tmp_path)}
    assert variantflow.parity(run) == {"available": False, "direct": True, "flow": False}
    (f / "variantcall.tsv").write_text(head + "g1\t20\tC\tT\t0\t1\t17\t1.0\ng1\t10\tA\tG\t1\t1\t17\t1.0\n")          # same calls, other order
    assert variantflow.parity(run)["same"] is True
    (f / "variantcall.tsv").write_text(head + "g1\t10\tA\tG\t1\t0\t\t\ng1\t30\tG\tA\t1\t0\t\t\n")
    p = variantflow.parity(run)
    assert p["same"] is False and p["only_direct"] == 1 and p["only_flow"] == 1 and p["differ"] == 1
    assert variantflow.available(run) and not variantflow.available({"run_dir": str(tmp_path / "x")})


def test_snakefile_collects_the_same_columns_as_the_direct_implementation():
    """The flow's scripts and variantcall.py must write the same table; a cheap structural check that the Snakefile and scripts stay in step (the real run is compared by variantflow.parity)."""
    from pathlib import Path
    root = Path(__file__).resolve().parent.parent / "flows" / "variantcheck"
    assert "genome\\tpos\\tref\\talt\\tsnippy\\tlofreq\\tdp\\taf" in (root / "scripts" / "collect.py").read_text(encoding="utf-8")
    src = (Path(__file__).resolve().parent.parent / "genomecurator" / "variantcall.py").read_text(encoding="utf-8")
    assert "genome\\tpos\\tref\\talt\\tsnippy\\tlofreq\\tdp\\taf" in src
    for rule in ("one_genome", "snippy", "lofreq", "calls", "collect"):
        assert f"rule {rule}:" in (root / "Snakefile").read_text(encoding="utf-8")


def test_flowlogic_is_pure_and_gives_the_report_and_review_bed(tmp_path):
    import ast
    from pathlib import Path
    from genomecurator import flowlogic
    src = (Path(flowlogic.__file__)).read_text(encoding="utf-8")
    imported = {n.names[0].name.split(".")[0] if isinstance(n, ast.Import) else (n.module or "").split(".")[0] for n in ast.walk(ast.parse(src)) if isinstance(n, (ast.Import, ast.ImportFrom))}
    assert imported <= {"csv", "collections"}, imported                 # the Snakemake environment has no other package
    v = tmp_path / "v.tsv"
    v.write_text("genome\tpos\tref\talt\tsnippy\tlofreq\tdp\taf\n" + "".join(f"g{i}\t50\tA\tG\t0\t1\t17\t1.0\n" for i in range(4)) + "g0\t60\tC\tT\t1\t1\t17\t1.0\n")
    m = tmp_path / "m.tsv"
    m.write_text("sample\tpos\tref\talt\ng0\t60\tC\tT\n")
    judged = {"g0", "g1", "g2", "g3"}
    callers = {"mutprofile": flowlogic.read_mutations(m, judged), **flowlogic.read_calls(v, judged)}
    rep = flowlogic.sites_report(callers)
    assert rep["sites_disagree"] == 1 and rep["sites"][0]["pos"] == 50 and rep["sites"][0]["lofreq"] == 4 and rep["sites"][0]["mutprofile"] == 0
    assert flowlogic.review_bed(rep["all_sites"], "ref", 3) == "ref\t49\t50\n" and flowlogic.review_bed(rep["all_sites"], "ref", 5) == ""
    assert flowlogic.read_calls(v, {"g0"}, lambda p: p < 55)["lofreq"]["g0"] == {(50, "A", "G")}      # window: the ITR side is left out


def test_callcache_shares_calls_between_runs_per_sequence_and_reference(tmp_path, monkeypatch):
    from genomecurator import callcache, config
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    k = callcache.key("KJ642613.1", "ACGT" * 10, "4.6.0", "2.1.5")
    assert k != callcache.key("KJ642613.1", "ACGT" * 10, "4.6.1", "2.1.5") and k != callcache.key("KJ642613.1", "ACGA" * 10, "4.6.0", "2.1.5")     # other reference or tool version = other cache
    run1 = tmp_path / "run1"
    (run1).mkdir()
    (run1 / "private_X.tsv").write_text("private:X\t10\tA\tG\t1\t1\t17\t1.0\n")
    (run1 / "A1.tsv").write_text("")                                                                  # a genome without any call is a result too
    recs = [("private:X", "acgtacgt"), ("A1", "TTTT")]
    assert callcache.store(recs, k, run1) == 2 and callcache.store(recs, k, run1) == 0
    run2 = tmp_path / "run2"
    n = callcache.prefill([("OTHER_NAME", "ACGTACGT"), ("B2", "TTTT"), ("NEW", "GGGG")], k, run2)       # same sequence under other ids; one new genome
    assert n == 2 and (run2 / "OTHER_NAME.tsv").read_text() == "OTHER_NAME\t10\tA\tG\t1\t1\t17\t1.0\n" and (run2 / "B2.tsv").read_text() == "" and not (run2 / "NEW.tsv").exists()
    assert callcache.prefill([("OTHER_NAME", "ACGTACGT")], k, run2) == 0                                 # never over a file that is there


def test_mutation_sites_for_the_timeline_and_the_shared_menu():
    import pytest
    from pathlib import Path
    from genomecurator import lineages, web
    assert lineages.parse_sites("7226G, C18850C,7226G") == [(7226, "G"), (18850, "C")]            # a leading reference base is accepted, repeats are dropped
    assert lineages.parse_sites(["a12t"]) == [(12, "T")]
    for bad in ("zzz", "7226", "7226N", "12-3A"):
        with pytest.raises(ValueError):
            lineages.parse_sites(bad)
    tdir = Path(web.__file__).parent / "templates"
    base = (tdir / "base.html").read_text(encoding="utf-8")
    assert "function gcMutMenu" in base.split("</head>")[0] and "function gcMenuClose" in base                    # defined in <head>: page scripts call it while the page is parsed
    assert "gcMutMenu(ev,RID" in (tdir / "heatmap.html").read_text(encoding="utf-8") and "gcMutMenu(ev,RID" in (tdir / "snipit.html").read_text(encoding="utf-8")
    assert "chipMenu" in (tdir / "novel.html").read_text(encoding="utf-8") and "sites=" in base


def test_coldspot_model_uses_the_apobec_context_and_corrects_for_the_number_of_regions():
    import random
    from genomecurator import outbreak
    rng = random.Random(3)
    L, itr = 12000, 1000
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    sites = list(range(itr + 1, L - itr + 1))
    # 500 events spread evenly over the interior, none inside region B (a cold region), twice as many inside region C (a hot one)
    B, C = (3001, 4000), (6001, 7000)
    events = []
    for p in rng.sample([x for x in sites if not (B[0] <= x <= B[1])], 500):
        events.append((seq[p - 1], p, "A" if seq[p - 1] != "A" else "C"))
    regions = [("A", 1001, 2000, [(1001, 2000)]), ("B", *B, [B]), ("C", *C, [C]), ("D", 8001, 9000, [(8001, 9000)])]
    rows = {r["name"]: r for r in outbreak.coldspot_model(events, seq, regions, itr)}
    assert rows["B"]["events"] == 0 and rows["B"]["expected"] > 5 and rows["B"]["coldspot"] is True and rows["B"]["p_bonferroni"] < 0.05
    assert not rows["A"]["coldspot"] and not rows["D"]["coldspot"] and rows["C"]["coldspot"] is False
    # too little expected: zero events is not surprising, so no coldspot however small p looks
    tiny = outbreak.coldspot_model(events[:3], seq, [("T", 3001, 3100, [(3001, 3100)])], itr)[0]
    assert tiny["expected"] < 5 and tiny["coldspot"] is False
    # the ITR is left out: a region wholly inside it is not tested
    assert outbreak.coldspot_model(events, seq, [("ITR", 1, 900, [(1, 900)])], itr) == []
    assert outbreak._pois_cdf(0, 2.0) == pytest_approx(0.1353352832366127) and outbreak._pois_cdf(-1, 2.0) == 0.0


def pytest_approx(x):
    import pytest
    return pytest.approx(x)


def test_the_page_says_when_the_code_is_newer_than_the_running_server(tmp_path, monkeypatch):
    import time
    from fastapi.testclient import TestClient
    from genomecurator import config, web
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    c = TestClient(web.create_app(scheduler=False))
    monkeypatch.setitem(web._stale_cache, "t", 0.0)
    monkeypatch.setattr(web, "_STARTED", time.time() + 3600)                       # started "in the future": nothing is newer
    assert "앱을 다시 시작하세요" not in c.get("/").text
    monkeypatch.setitem(web._stale_cache, "t", 0.0)
    monkeypatch.setattr(web, "_STARTED", 1.0)                                       # started long ago: every file is newer
    assert "앱을 다시 시작하세요" in c.get("/").text


def test_checkv_summary_is_read_into_qc_statuses_and_db_is_found(tmp_path, monkeypatch):
    from genomecurator import config, qc
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.delenv("GC_CHECKV_DB", raising=False)
    assert qc.checkv_db() is None
    (tmp_path / "checkv" / "checkv-db-v1.5" / "genome_db").mkdir(parents=True)
    assert qc.checkv_db().name == "checkv-db-v1.5"
    head = "contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n"
    f = tmp_path / "quality_summary.tsv"
    f.write_text(head + "A1\t197000\tNo\tNA\t190\t180\t0\tComplete\tHigh-quality\t100.0\tAAI-based\t0.0\t1.0\t\n"
                 "A2\t150000\tNo\tNA\t150\t140\t2\tLow-quality\tLow-quality\t61.5\tAAI-based\t3.2\t1.0\thigh kmer_freq may indicate large duplication\n"
                 "A3\t1000\tNo\tNA\t1\t0\t0\tNot-determined\tNot-determined\tNA\tNA\tNA\tNA\tno viral genes detected\n")
    r = qc.parse_checkv_summary(f)
    assert r["A1"]["status"] == "good" and r["A1"]["reasons"] == [] and r["A1"]["score"] == 100.0
    assert r["A2"]["status"] == "bad" and any("contamination 3.2%" in x for x in r["A2"]["reasons"]) and any("kmer_freq" in x for x in r["A2"]["reasons"])
    assert r["A3"]["status"] == "mediocre" and r["A3"]["metrics"]["checkv_quality"] == "Not-determined"


def test_dataset_made_from_a_selection_has_exactly_those_genomes(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient
    from genomecurator import datasets, private, sources, web
    _ds_env(tmp_path, monkeypatch)
    monkeypatch.setattr(sources, "pathoplexus_slug_for", lambda acc: None)
    private.import_texts(1, *_priv_files([{"PPX_accession": "LAB_V1", "strain": "kr-1", "date": "2025-02-03", "country": "Uganda", "clade": "Ib"},
                                          {"PPX_accession": "LAB_V2", "strain": "pe-1", "date": "2025", "country": "Peru"}], {"LAB_V1": "ACGT" * 40, "LAB_V2": "ACGT" * 41}))
    c = TestClient(web.create_app(scheduler=False))
    body = {"watch": 1, "name": "picked", "rules": {"only_manual": True, "require_sequence": True}, "manual_include": ["A1", "LAB_V2"], "include_private": True, "build": "clade-i"}
    ds = c.post("/api/datasets", json=body, headers={"origin": "http://testserver"}).json()
    assert "error" not in ds and ds["rules"]["only_manual"] is True and "고른 유전체만" in datasets.describe(ds["rules"])
    pub, priv = datasets.member_ids(1, datasets.get(ds["id"]))
    assert pub == {"A1"} and priv == {"LAB_V2"}                                  # nothing else matches although the rules name no clade or country
    nopriv = datasets.save(1, "picked public", {"only_manual": True}, manual_include=["A1", "LAB_V2"])               # private ones only when the dataset asks for them
    assert datasets.member_ids(1, nopriv) == ({"A1"}, set())
    assert datasets.member_ids(1, datasets.save(1, "all", {}))[0] >= {"A1"}      # the ordinary rules are unchanged
    page = c.get(f"/dataset/{ds['id']}")
    assert page.status_code == 200 and "dsAnalysis" in page.text and "makeDs" in c.get("/genomes?watch=1").text
    assert c.get(f"/api/datasets/{ds['id']}/nextstrain").status_code == 200


def test_pipeline_runs_the_steps_in_order_and_one_failure_does_not_stop_the_others(tmp_path, monkeypatch):
    from genomecurator import config, datasets, mutprofile, nextstrain, pipeline, squirrel, variantflow
    from genomecurator.db import connect
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    con = connect()
    con.execute("INSERT INTO watchlist(id,name,taxon) VALUES(1,'Orthopoxvirus monkeypox','Orthopoxvirus monkeypox')")
    con.execute("INSERT INTO watchlist(id,name,taxon) VALUES(2,'Henipavirus nipahense','Henipavirus nipahense')")
    con.commit(); con.close()
    calls = []
    monkeypatch.setattr(datasets, "run_qc", lambda ds, p, force: calls.append("qc") or {"public": {"sequences": 5}})
    monkeypatch.setattr(squirrel, "run", lambda ds, cores, force, note, p: calls.append("squirrel") or {"id": 10, "reused": True})
    monkeypatch.setattr(mutprofile, "run", lambda ds, cores, force, note, p: calls.append("mutprofile") or (_ for _ in ()).throw(RuntimeError("alignment failed")))
    monkeypatch.setattr(nextstrain, "run", lambda ds, rt, cores, force, note, p: calls.append("nextstrain") or {"id": 12})
    monkeypatch.setattr(variantflow, "available", lambda r: False)
    monkeypatch.setattr(variantflow, "run_recorded", lambda r, cores, p: calls.append("variants") or {"run": 99, "genomes": 3, "from_cache": 1})
    prog = []
    out = pipeline.run_all({"id": 1, "watch_id": 1, "name": "d"}, None, False, None, lambda pct, m="": prog.append(pct))
    assert calls == ["qc", "squirrel", "mutprofile", "nextstrain"]                                       # mutprofile failed: variants are skipped, nextstrain still runs
    s = out["steps"]
    assert s["qc"]["status"] == "done" and s["squirrel"]["status"] == "reused" and s["mutprofile"]["status"] == "error" and "alignment failed" in s["mutprofile"]["detail"]
    assert s["variants"]["status"] == "skipped" and s["nextstrain"]["status"] == "done"
    assert prog[-1] == 100 and prog == sorted(prog)
    # a virus without analyses: QC only, and the user is told
    calls.clear()
    out2 = pipeline.run_all({"id": 2, "watch_id": 2, "name": "n"}, None, False, None)
    assert calls == ["qc"] and list(out2["steps"]) == ["qc"] and any("QC만" in n for n in out2["notes"])
    # chosen steps only
    calls.clear()
    pipeline.run_all({"id": 1, "watch_id": 1, "name": "d"}, ["qc", "nextstrain"])
    assert calls == ["qc", "nextstrain"]


def test_qc_headline_follows_the_tool_priority():
    from genomecurator import qc
    assert qc.PRIORITY[:3] == ("squirrel", "nextclade", "checkv")
    assert qc.overall({"squirrel": "good", "nextclade": "bad", "checkv": "bad", "seqstats": "bad"}) == "good"          # Squirrel counts first
    assert qc.overall({"nextclade": "mediocre", "checkv": "bad"}) == "mediocre"                                           # then Nextclade
    assert qc.overall({"checkv": "bad", "seqstats": "good"}) == "bad"                                                      # CheckV before the built-in statistics
    assert qc.overall({"squirrel": "error", "nextclade": "good"}) == "good"                                               # a tool that failed has no verdict
    assert qc.overall({"squirrel": "error"}) == "error" and qc.overall({}) is None


def test_clock_inputs_outputs_and_group_tmrca(tmp_path):
    import pytest
    from genomecurator import clock
    s = clock.parse_settings({"source": "squirrel", "clock_filter": "2.5", "min_group": 3})
    assert s["clock_filter"] == 2.5 and s["min_group"] == 3 and s["reroot"] == "least-squares" and s["covariation"] is True
    for bad in ({"source": "x"}, {"reroot": "x"}, {"clock_filter": 99}, {"coalescent": "foo"}, {"min_group": 1}, {"min_precision": "week"}):
        with pytest.raises(ValueError):
            clock.parse_settings(bad)
    assert clock.parse_settings({"coalescent": "skyline"})["coalescent"] == "skyline" and clock.parse_settings({"coalescent": "0.5"})["coalescent"] == "0.5"
    # dates: a day as it is, a month as its range, nothing below the chosen precision
    assert clock._treetime_date("2025-03-04", "day") == ("2025-03-04", "day")
    assert clock._treetime_date("2025-02", "day") is None and clock._treetime_date("2025-02", "month") == ("[2025-02-01:2025-02-28]", "month")
    assert clock._treetime_date("2024-12", "month")[0] == "[2024-12-01:2024-12-31]" and clock._treetime_date("2025", "month") is None
    assert clock._treetime_date("2025", "year") == ("[2025-01-01:2025-12-31]", "year") and clock._treetime_date("", "year") is None
    assert abs(clock.given_decimal("2025-07-02") - 2025.5) < 0.01 and abs(clock.given_decimal("[2025-01-01:2025-12-31]") - 2025.5) < 0.01
    # TreeTime's files
    (tmp_path / "mc.txt").write_text("Root-Tip-Regression:\n --rate:\t2.157e-05 +/- 3.57e-06 (one std-dev)\n --r^2:  \t0.33\n")
    assert clock.parse_clock_txt(tmp_path / "mc.txt") == {"rate": 2.157e-05, "rate_sd": 3.57e-06, "r2": 0.33}
    (tmp_path / "rtt.csv").write_text("#Dates of nodes ...\nname, date, root-to-tip distance, clock-deviation\nA, 2025.1, 0.0001, -0.00002\nNODE_0000001, 2024.0, 0.00005, 0.0\nB, 2025.9, 0.0003, 0.0001\n")
    (tmp_path / "out.tsv").write_text("\tgiven_date\tapparent_date\tresidual\nB\t2025.9\t2031.0\t5.1\n")
    pts, outl = clock.parse_rtt(tmp_path / "rtt.csv", tmp_path / "out.tsv")
    assert [p["name"] for p in pts] == ["A", "B"] and outl == {"B"}
    # the MRCA of a group in a time tree: a monophyletic group, one that is not, and one that is too small
    def node(name, date, lo, hi, children=None):
        n = {"name": name, "node_attrs": {"num_date": {"value": date, "confidence": [lo, hi]}}}
        if children:
            n["children"] = children
        return n
    tree = node("N0", 2020.0, 2019.0, 2021.0, [node("N1", 2023.0, 2022.5, 2023.4, [node("a1", 2024.0, 2024.0, 2024.0), node("a2", 2024.5, 2024.5, 2024.5)]),
                                               node("N2", 2022.0, 2021.5, 2022.5, [node("b1", 2025.0, 2025.0, 2025.0), node("c1", 2025.2, 2025.2, 2025.2)])])
    g = {x["group"]: x for x in clock.group_tmrca(tree, {"A": {"a1", "a2"}, "BC": {"a1", "b1"}, "small": {"c1"}}, 2)}
    assert g["A"]["tmrca"] == 2023.0 and g["A"]["extra"] == 0 and g["A"]["lo"] == 2022.5
    assert g["BC"]["tmrca"] == 2020.0 and g["BC"]["extra"] == 2 and "small" not in g                       # not monophyletic: the root, with the other genomes counted
    w = clock.warnings({"clock": {"r2": 0.1}, "span_years": 2.9, "time": {"root": {"lo": 2019.0, "hi": 2022.4}}, "n_outliers": 20, "n_dated": 100, "precision": {"month": 50}, "n_dropped_no_date": 3})
    assert len(w) >= 6 and any("R² 0.10" in x for x in w) and any("2.9년" in x for x in w) and any("외삽" in x for x in w)


def test_beast_package_subsample_xml_and_posterior_reading(tmp_path):
    import random
    import xml.etree.ElementTree as ET
    import pytest
    from genomecurator import beastprep
    s = beastprep.parse_settings({"n_tips": 20, "chain_length": 50000, "samples": 500})
    assert s["clock"] == "relaxed" and s["tree_prior"] == "exponential"
    for bad in ({"n_tips": 5}, {"clock": "x"}, {"tree_prior": "x"}, {"samples": 5}, {"chain_length": 10}):
        with pytest.raises(ValueError):
            beastprep.parse_settings(bad)
    rng = random.Random(2)
    tips = [{"name": f"g{i}", "date": 2023.0 + (i % 36) / 12 + rng.random() / 12, "outlier": i % 17 == 0, "precision": "range" if i % 9 == 0 else "day", "lineage": f"A.{i % 4}", "clade": "Ib", "country": "X"} for i in range(200)]
    ch = beastprep.choose_tips(tips, 30, seed=1)
    assert len(ch) == 30 and all(not t["outlier"] and t["precision"] == "day" for t in ch)
    assert ch == beastprep.choose_tips(tips, 30, seed=1) and [t["date"] for t in ch] == sorted(t["date"] for t in ch)                  # deterministic, in date order
    assert max(t["date"] for t in ch) - min(t["date"] for t in ch) > 2.5 and {t["lineage"] for t in ch} == {"A.0", "A.1", "A.2", "A.3"}   # the whole period and every lineage kept
    assert len(beastprep.choose_tips(tips[:12], 30)) <= 12
    seqs = {t["name"]: "".join(rng.choice("ACGT") for _ in range(60)) for t in ch}
    for clk in ("strict", "relaxed"):
        for tp in ("constant", "exponential"):
            xml = beastprep.build_xml(seqs, {t["name"]: t["date"] for t in ch}, {"lineage.A.1": [t["name"] for t in ch if t["lineage"] == "A.1"]}, {**s, "clock": clk, "tree_prior": tp}, 2e-5)
            root = ET.fromstring(xml)                                                                                                    # well-formed
            ids = {e.get("id") for e in root.iter() if e.get("id")}
            assert {"taxa", "alignment", "treeModel", "age(root)", "age(lineage.A.1)", "mcmc", "fileLog", "treeLikelihood", "demographic"} <= ids
            assert ("ucld.mean" in ids) == (clk == "relaxed") and ("clock.rate" in ids) == (clk == "strict") and ("demographic.growthRate" in ids) == (tp == "exponential")
            assert len(list(root.iter("taxon"))) > 2 * len(ch) and root.find("mcmc").get("chainLength") == "50000"
    # the posterior: a log with a known answer
    log = tmp_path / "b.log"
    vals = [2023.0 + 0.5 * math_sin(i) for i in range(400)]
    log.write_text("# comment\nstate\tposterior\tage(root)\tclock.rate\n" + "".join(f"{i * 100}\t-1.0\t{v:.5f}\t{1e-5 + i * 1e-9:.3e}\n" for i, v in enumerate(vals)))
    r = beastprep.read_log(log)
    a = r["params"]["age(root)"]
    assert r["states"] == 400 and r["last_state"] == 39900 and 2022.5 <= a["lo"] <= a["median"] <= a["hi"] <= 2023.5 and a["ess"] >= 1
    lo, hi = beastprep.hpd(list(range(100)), 0.95)
    assert (lo, hi) == (0, 94) or (hi - lo) == 94
    assert beastprep.ess([float(i % 2) for i in range(400)]) > 100 and beastprep.ess([float(i) for i in range(400)]) < 20     # independent-looking vs a drifting chain


def math_sin(i):
    import math
    return math.sin(i * 0.9)


def test_bets_xml_settings_and_verdict():
    import random
    import xml.etree.ElementTree as ET
    import pytest
    from genomecurator import beastprep
    for bad in ({"n_tips": 5}, {"path_steps": 2}, {"step_chain": 10}, {"initial_chain": 100}, {"clock": "x"}):
        with pytest.raises(ValueError):
            beastprep.parse_bets_settings(bad)
    s = beastprep.parse_bets_settings({})
    assert s["path_steps"] >= 50 and s["step_chain"] >= 250000                                              # the defaults follow the recommendation for a reliable estimate
    rng = random.Random(4)
    names = [f"s{i}" for i in range(12)]
    seqs = {n: "".join(rng.choice("ACGT") for _ in range(40)) for n in names}
    dates = {n: 2024 + i / 10 for i, n in enumerate(names)}
    ms = {"n_tips": 12, "chain_length": 30000, "samples": 300, "clock": "strict", "tree_prior": "exponential", "seed": 1, "threads": 0}
    xml = beastprep.build_xml(seqs, dates, {}, ms, 1e-5, mle={"chain_length": 5000, "path_steps": 20, "alpha": 0.3}, proper=True)
    root = ET.fromstring(xml)
    mle = root.find("marginalLikelihoodEstimator")
    assert mle.get("pathSteps") == "20" and mle.get("chainLength") == "5000" and mle.get("pathScheme") == "betaquantile"
    assert root.find("pathSamplingAnalysis") is not None and root.find("steppingStoneSamplingAnalysis") is not None
    assert "oneOnXPrior" not in xml and 'gammaPrior shape="0.001"' not in xml                                       # proper priors only: a marginal likelihood needs them
    vague = beastprep.build_xml(seqs, dates, {}, ms, 1e-5)
    assert "oneOnXPrior" in vague and "marginalLikelihoodEstimator" not in vague                                     # the ordinary analysis is unchanged
    same = beastprep.build_xml(seqs, {n: 2026.0 for n in names}, {}, ms, 1e-5, mle={"chain_length": 5000, "path_steps": 20}, proper=True)
    assert len({e.find("date").get("value") for e in ET.fromstring(same).find("taxa")}) == 1                          # the isochronous analysis: one date for every genome
    assert [beastprep.bets_verdict(x)[0] for x in (12.0, 4.0, 2.0, 0.3, -3.0, None)] == ["strong", "positive", "weak", "none", "against", "unknown"]
    log_line = "log marginal likelihood (using path sampling) from pathLikelihood.delta = -5003.54\nlog marginal likelihood (using stepping stone sampling) from pathLikelihood.delta = -5004.2\n"
    assert float(beastprep._MLE_RE["ps"].search(log_line)[1]) == -5003.54 and float(beastprep._MLE_RE["ss"].search(log_line)[1]) == -5004.2


def test_gss_xml_has_a_working_prior_for_every_free_parameter():
    import random
    import xml.etree.ElementTree as ET
    from genomecurator import beastprep
    rng = random.Random(6)
    names = [f"s{i}" for i in range(10)]
    seqs = {n: "".join(rng.choice("ACGT") for _ in range(50)) for n in names}
    dates = {n: 2024 + i / 10 for i, n in enumerate(names)}
    for clock_, prior in (("strict", "constant"), ("relaxed", "exponential")):
        ms = {"n_tips": 10, "chain_length": 2000, "samples": 100, "clock": clock_, "tree_prior": prior, "seed": 1, "threads": 0}
        xml = beastprep.build_xml(seqs, dates, {}, ms, 1e-5, "gss", proper=True, fixed_freq=True, init={"kappa": 3.2, "clock.rate": 2e-5, "ucld.mean": 3e-5},
                                  gss={"log": "beast.log", "burnin": 500, "chain_length": 1000, "path_steps": 10, "alpha": 0.3})
        root = ET.fromstring(xml)
        wp = root.find("marginalLikelihoodEstimator").find("pathLikelihood").find("destination").find("workingPrior")
        covered = {e.find("parameter").get("idref") for e in wp if e.find("parameter") is not None}
        free = {"kappa", "alpha", "demographic.popSize"} | ({"clock.rate"} if clock_ == "strict" else {"ucld.mean", "ucld.stdev"}) | ({"demographic.growthRate"} if prior == "exponential" else set())
        assert covered == free, (clock_, prior, covered ^ free)                                                  # a free parameter without a working prior would break the estimate
        assert wp.find("coalescentLikelihood").get("idref") == "coalescentReference" and root.find("generalizedSteppingStoneSamplingAnalysis") is not None
        assert all(e.get("fileName") == "beast.log" and e.get("burnin") == "500" for e in wp if e.get("parameterColumn"))
        assert "dirichletPrior" not in xml and "deltaExchange" not in xml                                         # fixed base frequencies: no frequency simplex to give a working prior
        freq = [e for e in root.iter("parameter") if e.get("id") == "frequencies"][0].get("value").split()
        assert len(freq) == 4 and abs(sum(map(float, freq)) - 1) < 1e-3 and freq != ["0.25"] * 4                 # the alignment's own
        assert ('id="kappa" value="3.2"' in xml) and (('id="clock.rate" value="2e-05"' in xml) if clock_ == "strict" else ('id="ucld.mean" value="3e-05"' in xml))   # the chain starts inside the working prior
    assert "workingPrior" not in beastprep.build_xml(seqs, dates, {}, ms, 1e-5, mle={"chain_length": 1000, "path_steps": 5}, proper=True)
    s = beastprep.parse_bets_settings({})
    assert s["estimator"] == "both" and beastprep.parse_bets_settings({"estimator": "gss"})["estimator"] == "gss"
    import pytest
    with pytest.raises(ValueError):
        beastprep.parse_bets_settings({"estimator": "harmonic"})
    assert float(beastprep._MLE_RE["gss"].search("log marginal likelihood (using generalized stepping stone sampling) from (pathLikelihood.source - pathLikelihood.destination) = -1220.4")[1]) == -1220.4


def test_mutation_gantt_orders_by_first_observation_and_counts_per_period(tmp_path, monkeypatch):
    import csv
    import pytest
    from genomecurator import insight, mutgantt, mutprofile
    # three mutations carried by genomes collected on known dates; M2 first observed earliest, M3 only from a genome with a month-only date
    meta = {"g1": ("2025-01-05", "A"), "g2": ("2025-01-20", "A"), "g3": ("2025-03-02", "B"), "g4": ("2025-03", "B"), "g5": ("", "B")}
    monkeypatch.setattr(mutprofile, "_run_meta", lambda run: {k: {"clade": "Ib", "lineage": l, "country": "X", "date": d, "outbreak": ""} for k, (d, l) in meta.items()})
    monkeypatch.setattr(mutprofile, "summary", lambda run: {"reference": ""})
    monkeypatch.setattr(mutprofile, "reference_by_name", lambda ref: (None, None))
    monkeypatch.setattr(mutprofile, "read_fasta", lambda fa: ("", ""))
    def row(s, pos, alt, apo=""):
        return {"sample": s, "pos": str(pos), "ref": "A", "alt": alt, "gene": "G1" if pos < 200 else "", "aa_change": "", "effect": "synonymous", "profile_of_interest": apo}
    rows = [row("g3", 100, "T"), row("g4", 100, "T"), row("g5", 100, "T"),                       # M1 (100T): first seen 2025-03-02
            row("g1", 300, "C", "apobec"), row("g2", 300, "C", "apobec"), row("g3", 300, "C", "apobec"),  # M2 (300C): first seen 2025-01-05
            row("g4", 50, "G"), row("g5", 50, "G"), row("g3", 50, "G")]                                 # M3 (50G): g4 has only a month
    monkeypatch.setattr(insight, "_mutations", lambda run, max_n=None: rows)
    g = mutgantt.build({}, "month", min_carriers=2)
    assert [x["key"] for x in g["rows"]] == ["300C", "50G", "100T"] and g["total_rows"] == 3                      # earliest first observation on top
    by = {x["key"]: x for x in g["rows"]}
    assert by["300C"]["first"] == "2025-01-05" and by["300C"]["last"] == "2025-03-02" and by["300C"]["cells"] == {"2025-01": 2, "2025-03": 1} and by["300C"]["apobec"] == 1
    assert by["100T"]["n"] == 3 and by["100T"]["n_dated"] == 2 and by["100T"]["first"] == "2025-03"               # g5 has no date: counted as a carrier, not placed in time; a month-only date counts as the start of its month, so it is the earliest possible
    assert g["periods"] == ["2025-01", "2025-02", "2025-03"] and g["no_date"] == 1 and g["month_only"] == 1
    assert by["300C"]["by_lineage"] == {"A": 2, "B": 1} and by["300C"]["by_country"] == {"X": 3} and by["300C"]["by_continent"] == {}      # carriers per lineage / country / continent, for the highlight legend
    w = mutgantt.build({}, "week", min_carriers=2)
    assert w["periods"][0] == "2025-W01" and w["rows"][0]["key"] == "300C" and sum(w["rows"][0]["cells"].values()) == 3
    assert {sum(x["cells"].values()) for x in w["rows"] if x["key"] == "50G"} == {1}                              # the month-only genome is not on a week
    assert [x["key"] for x in mutgantt.build({}, "month", min_carriers=2, sort="pos")["rows"]] == ["50G", "100T", "300C"]
    assert [x["key"] for x in mutgantt.build({}, "month", min_carriers=2, sort="n")["rows"]] == ["300C", "50G", "100T"]               # equal counts: the earlier first observation first
    assert [x["key"] for x in mutgantt.build({}, "month", min_carriers=2, apobec_only=True)["rows"]] == ["300C"]
    assert [x["key"] for x in mutgantt.build({}, "month", min_carriers=2, start=1, end=120)["rows"]] == ["50G", "100T"]
    assert mutgantt.build({}, "month", min_carriers=4)["total_rows"] == 0                                         # a mutation needs enough carriers
    assert [x["key"] for x in mutgantt.build({}, "month", min_carriers=2, lineage_set=["A"])["rows"]] == ["300C"]
    assert len(mutgantt.build({}, "month", min_carriers=2, max_rows=2)["rows"]) == 2
    with pytest.raises(ValueError):
        mutgantt.build({}, "month", gene="NOPE")
    # lineage-defining mutations: tagged, and filtered with defining=only/none
    monkeypatch.setattr(mutgantt, "defining_map", lambda run, meta, fa: {(300, "C"): ["A"]})
    d = mutgantt.build({}, "month", min_carriers=2)
    assert {x["key"]: x["defines"] for x in d["rows"]} == {"300C": ["A"], "50G": [], "100T": []} and d["n_defining_shown"] == 1 and d["definitions"]
    assert [x["key"] for x in mutgantt.build({}, "month", min_carriers=2, defining="only")["rows"]] == ["300C"]
    assert [x["key"] for x in mutgantt.build({}, "month", min_carriers=2, defining="none")["rows"]] == ["50G", "100T"]
    # one mutation over time: carriers per period by lineage against the genomes sampled in that period
    sr = mutgantt.series({}, 300, "C", "month")
    assert sr["periods"] == ["2025-01", "2025-02", "2025-03"] and sr["sampled"] == [2, 0, 2] and sr["carriers"] == [2, 0, 1] and sr["n_carriers"] == 3 and sr["undated"] == 1
    assert sr["lineages"][0]["lineage"] == "A" and sr["lineages"][0]["counts"] == [2, 0, 0] and sr["lineages"][1]["lineage"] == "B" and sr["lineages"][1]["counts"] == [0, 0, 1]
    sw = mutgantt.series({}, 100, "T", "week")
    assert sum(sw["carriers"]) == 1 and len(sw["periods"]) == 1 and sum(sw["sampled"]) == 1                                                 # the month-only genome is left out of the weekly view
    assert mutgantt.series({}, 999, "A", "month")["n_carriers"] == 0


def test_empirical_base_frequencies_sum_to_one_for_beast():
    import random
    import xml.etree.ElementTree as ET
    from genomecurator import beastprep
    for seed in range(40):                                                   # rounding to few decimals made some sets sum to 1.000001, which BEAST refuses
        rng = random.Random(seed)
        w = [rng.random() + 0.1 for _ in range(4)]
        seqs = {f"s{i}": "".join(rng.choices("ACGT", weights=w, k=777)) for i in range(5)}
        xml = beastprep.build_xml(seqs, {n: 2024.0 for n in seqs}, {}, {"n_tips": 5, "chain_length": 1000, "samples": 100, "clock": "strict", "tree_prior": "constant", "seed": 1, "threads": 0}, 1e-5, fixed_freq=True)
        f = [float(x) for x in [e for e in ET.fromstring(xml).iter("parameter") if e.get("id") == "frequencies"][0].get("value").split()]
        assert len(f) == 4 and abs(sum(f) - 1.0) < 1e-9


def test_adding_a_virus_that_is_already_watched_is_not_an_error(tmp_path, monkeypatch):
    import pytest
    from genomecurator import config, curate, sources
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    monkeypatch.setattr(sources, "ncbi_resolve_virus", lambda q: (3052458, "Orthoebolavirus bundibugyoense"))
    monkeypatch.setattr(sources, "ncbi_species_scope", lambda tid: None)
    first = curate.add_watch("Orthoebolavirus bundibugyoense")
    with pytest.raises(curate.AlreadyWatched) as e:                              # same name: caught before the database refuses it with a traceback
        curate.add_watch("Orthoebolavirus bundibugyoense")
    assert e.value.watch_id == first and e.value.name == "Orthoebolavirus bundibugyoense" and isinstance(e.value, ValueError)
    with pytest.raises(curate.AlreadyWatched) as e2:                             # another name for the same taxon
        curate.add_watch("Bundibugyo virus", "Orthoebolavirus bundibugyoense")
    assert e2.value.watch_id == first


def test_nextclade_calls_are_moved_to_the_run_reference_and_compared(tmp_path, monkeypatch):
    from genomecurator import nextcladecall, variantcall
    # dataset reference (NC) and run reference differ at NC 5 (T) / run 5 (A) and the run reference has one extra base up front (run position = NC position + 1)
    nc = "ACGTTACGTA"
    run = "G" + "ACGTAACGTA"                                  # run pos 6 = NC pos 5: run A vs NC T; everything else equal
    to_nc = {p: p - 1 for p in range(2, 12)}
    tsv = tmp_path / "nc.tsv"
    tsv.write_text("seqName\tsubstitutions\tdeletions\tmissing\tnonACGTs\terrors\n"
                   "g1\tG3T\t\t\t\t\n"                       # NC 3 G->T: run pos 4 is G -> call (4, G, T); NC 5 not listed -> the genome has T, the run reference A -> call (6, A, T)
                   "g2\t\t\t5\t\t\n"                          # NC 5 missing -> unknown, no call there
                   "g3\t\t\t\t\tbad sequence\n", encoding="utf-8")
    calls = nextcladecall.parse_calls(tsv, nc, run, to_nc)
    assert calls["g1"] == {(4, "G", "T"), (6, "A", "T")} and calls["g2"] == set() and "g3" not in calls
    assert nextcladecall.parse_calls(tsv, nc, run, to_nc, window=lambda p: p != 6)["g1"] == {(4, "G", "T")}      # the ITR window is honoured
    assert nextcladecall._ranges("N:5-9,R:12") == [(5, 9), (12, 12)] and nextcladecall._ranges("") == []
    # one mutation tool by tool
    callers = {"mutprofile": {"g1": {(4, "G", "T")}, "g2": {(4, "G", "T")}, "g3": set()},
               "snippy": {"g1": {(4, "G", "T")}, "g2": set(), "g3": set()},
               "nextclade": {"g1": {(4, "G", "T")}, "g2": {(4, "G", "T")}}}
    monkeypatch.setattr(variantcall, "_callers", lambda run, max_n, source="direct": callers)
    sc = variantcall.site_calls({}, 4, "T")
    assert sc["per_caller"]["mutprofile"] == {"judged": 3, "called": 2} and sc["per_caller"]["nextclade"] == {"judged": 2, "called": 2} and sc["ref"] == "G"
    pats = {tuple(r["by"]): r for r in sc["patterns"]}
    assert pats[("mutprofile", "snippy", "nextclade")]["all"] and pats[("mutprofile", "snippy", "nextclade")]["genomes"] == ["g1"]
    assert not pats[("mutprofile", "nextclade")]["all"] and pats[("mutprofile", "nextclade")]["missed_by"] == ["snippy"]       # g2: snippy did not call it
    assert sc["n_called"] == 2 and sc["n_agree"] == 1 and sc["tool_specific"] == []


def test_lineage_tree_and_growth(monkeypatch):
    import numpy as np
    from genomecurator import lineages, lineagetree
    defs = {"sh1": {"": {"name": "", "parent": "", "clade": "Ib", "designation_date": "2023-01-01", "snps": [(10, "A")]},
                    "A": {"name": "A", "parent": "None", "clade": "Ib", "designation_date": "2023-02-01", "snps": [(10, "A"), (20, "C")]},
                    "A.1": {"name": "A.1", "parent": "A", "clade": "Ib", "designation_date": "2023-03-01", "snps": [(10, "A"), (20, "C"), (30, "G")]},
                    "A.2": {"name": "A.2", "parent": "A", "clade": "Ib", "designation_date": "2023-03-01", "snps": [(10, "A"), (20, "C"), (40, "T")]}}}
    monkeypatch.setattr(lineages, "definitions", lambda: defs)
    monkeypatch.setattr(lineages, "commit", lambda: "abc")
    monkeypatch.setattr(lineages, "nc_reference", lambda: ("NC", "ACGT" * 20, []))
    # a clean planted series: per ISO week w (starting 2025-01-06) A.2 has 30 - 2w genomes, A.1 has 4 + 2w
    gs = []
    from datetime import date, timedelta
    for w in range(10):
        d0 = date(2025, 1, 6) + timedelta(days=7 * w)
        gs += [{"id": f"x{w}-{i}", "lineage": "A.2", "outbreak": "sh1", "clade": "Ib", "country": "X", "date": d0.isoformat()} for i in range(30 - 2 * w)]
        gs += [{"id": f"y{w}-{i}", "lineage": "A.1", "outbreak": "sh1", "clade": "Ib", "country": "Y", "date": d0.isoformat()} for i in range(4 + 2 * w)]
    gs += [{"id": "r", "lineage": "", "outbreak": "sh1", "clade": "Ib", "country": "Y", "date": ""}]
    monkeypatch.setattr(lineages, "_genomes", lambda run: gs)
    t = lineagetree.tree({})
    nd = {n["id"]: n for n in t["nodes"]}
    assert nd["sh1/A"]["parent"] == "ob:sh1" and nd["sh1/A.1"]["parent"] == "sh1/A" and nd["ob:sh1"]["parent"] == "Ib" and nd["Ib"]["parent"] == "I"      # "None" as a parent means the outbreak root
    assert [s["label"] for s in nd["sh1/A.1"]["own_snps"]] == ["A4G"] or len(nd["sh1/A.1"]["own_snps"]) == 1
    assert nd["sh1/A.1"]["cum_snps"] == 3 and nd["sh1/A"]["cum_snps"] == 2 and nd["ob:sh1"]["cum_snps"] == 1                                        # mutations accumulate down the path
    assert nd["mpox"]["n_sub"] == len(gs) and nd["sh1/A.1"]["n_own"] == sum(1 for g in gs if g["lineage"] == "A.1")
    g = lineagetree.growth({}, "week", 20, min_n=10)
    r = {x["lineage"]: x for x in g["rows"]}
    assert g["reference"] in ("A.1", "A.2") and set(r) == {"A.1", "A.2"} and g["genomes"] == sum(1 for x in gs if x["date"] and x["lineage"])
    up = r["A.1"]["vs_ref"] if g["reference"] == "A.2" else -r["A.2"]["vs_ref"]
    assert up > 0.1 and abs(sum(x["share_last"] for x in g["rows"]) - 1) < 1e-9                                                                    # A.1 gains on A.2, the fitted shares sum to one
    lo, hi = (r["A.1"]["vs_ref_lo"], r["A.1"]["vs_ref_hi"]) if g["reference"] == "A.2" else (-r["A.2"]["vs_ref_hi"], -r["A.2"]["vs_ref_lo"])
    assert lo > 0 and hi > lo                                                                                                                  # a clear trend: the interval excludes zero
    assert "error" in lineagetree.growth({}, "week", 20, country="Z")                                                                           # no genomes for that country
    only = lineagetree.growth({}, "week", 20, country="Y", min_n=10)
    assert "error" in only                                                                                                                      # one lineage alone cannot be compared


def test_geo_resolve_gives_the_continent():
    from genomecurator import geo
    assert geo.resolve("Uganda")["continent"] == "Africa" and geo.resolve("Germany")["continent"] == "Europe" and geo.resolve("USA")["continent"] == "North America"
    assert geo.resolve("Singapore")["continent"] == "Asia" and geo.resolve("Réunion")["continent"] == "Africa"                 # places the 110m outlines leave out
    assert all(geo.POINT_CONTINENT[v[0]] for v in geo.POINTS.values())


def test_continent_comparison_separates_lineage_effects_from_real_differences(monkeypatch):
    import math
    from genomecurator import continents, insight, mutprofile
    assert abs(continents.chi2_sf(3.841458820694124, 1) - 0.05) < 1e-6 and abs(continents.chi2_sf(9.210340371976184, 2) - 0.01) < 1e-6
    q = continents.bh([0.001, 0.04, 0.03, 0.5])
    assert [round(x, 4) for x in q] == [0.004, 0.0533, 0.0533, 0.5]
    # 3 continents x 2 lineages: mutation 100T marks lineage L1 (which is mostly in Europe: a lineage effect), mutation 200G is carried by Africa only inside both lineages (a real difference)
    meta, rows = {}, []
    plan = [("Germany", "L1", 24, 6), ("Germany", "L2", 6, 12), ("Uganda", "L1", 6, 6), ("Uganda", "L2", 24, 12), ("Thailand", "L1", 6, 6), ("Thailand", "L2", 24, 12)]
    i = 0
    for country, lin, n, _x in plan:
        for k in range(n):
            s = f"g{i}"
            i += 1
            meta[s] = {"clade": "Ib", "lineage": lin, "country": country, "date": "2025-01-01", "outbreak": ""}
            if lin == "L1":
                rows.append({"sample": s, "pos": "100", "ref": "A", "alt": "T", "gene": "G1", "aa_change": "", "effect": "", "profile_of_interest": ""})
            if country == "Uganda" and k % 2 == 0:
                rows.append({"sample": s, "pos": "200", "ref": "C", "alt": "G", "gene": "G2", "aa_change": "", "effect": "", "profile_of_interest": ""})
    monkeypatch.setattr(mutprofile, "_run_meta", lambda run: meta)
    monkeypatch.setattr(mutprofile, "good_samples", lambda run, mx: None)
    monkeypatch.setattr(insight, "_mutations", lambda run, mx=None: rows)
    d = continents.compare({}, "country", min_group=5, min_carriers=5, adjust="lineage", perms=300)
    by = {x["key"]: x for x in d["rows"]}
    assert set(d["groups"]) == {"Germany", "Uganda", "Thailand"} and d["tested"] == 2
    assert by["100T"]["q"] < 0.05 and by["100T"]["q_adj"] > 0.05            # the difference is the lineage mix: gone once lineages are held fixed
    assert by["200G"]["q"] < 0.05 and by["200G"]["q_adj"] < 0.05            # Uganda only, in both lineages: stays
    assert by["200G"]["top"] == "Uganda" and d["significant_raw"] == 2 and d["significant_adjusted"] == 1
    assert "error" in continents.compare({}, "country", min_group=500)
    assert continents.group_of("USA", "continent") == "North America" and continents.group_of("", "continent") == ""


def test_intake_history_counts_intake_and_collection_dates(tmp_path, monkeypatch):
    from genomecurator import config, db, intake
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    db.init_db() if hasattr(db, "init_db") else None
    with db.connect() as con:
        con.execute("INSERT INTO watchlist(id,name,taxon) VALUES(1,'V','v')")
        rows = [("a", "2026-10-16T17:00:00", "2025-03-02"), ("b", "2026-10-16T17:00:00", "2025-03"), ("c", "2026-10-19T10:00:00", "2025-04-10"), ("d", "2026-10-19T10:00:00", "1958"), ("e", "2026-10-19T10:00:00", "")]
        for acc, ins, col in rows:
            con.execute("INSERT INTO genome(genome_id,watch_id,accession,date_inserted,collection_date) VALUES(?,?,?,?,?)", (acc, 1, acc, ins, col))
        con.execute("INSERT INTO changeset(id,watch_id,status,applied_at) VALUES(7,1,'APPLIED','2026-10-19T10:00:00')")
        con.execute("INSERT INTO genome_version(genome_id,snapshot_json,changeset_id,changed_at) VALUES('a','{}',7,'2026-10-19 10:00:01')")
    d = intake.history(1, "day", "month")
    assert d["periods"] == ["2026-10-16", "2026-10-17", "2026-10-18", "2026-10-19"] and d["added"] == [2, 0, 0, 3] and d["changed"] == [0, 0, 0, 1]
    assert d["applied"][0]["id"] == 7 and d["applied"][0]["period"] == "2026-10-19" and d["total"] == 5
    assert d["cperiods"] == ["2025-03", "2025-04"] and d["collected"] == [2, 1] and d["year_only"] == 1 and d["no_date"] == 1
    w = intake.history(1, "week", "week")
    assert w["added"] == [2, 3] and w["coarse_skipped"] == 1                        # the month-only date cannot be put on a week
    s = intake.history(1, "day", "month", sel="2026-10-19")
    assert s["sel_n"] == 3 and sum(s["collected"]) == 1 and s["year_only"] == 1


def test_new_genomes_for_an_existing_dataset(tmp_path, monkeypatch):
    from genomecurator import config, datasets, db, export, private
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    with db.connect() as con:
        con.execute("INSERT INTO watchlist(id,name,taxon) VALUES(1,'V','v')")
    ds = datasets.save(1, "Ib", {"clade": ["Ib"], "country": ["Uganda", "Rwanda"], "year_min": 2023, "require_sequence": True})
    base = datasets.baseline(ds)
    assert base[1] == "데이터셋 생성"
    g = lambda acc, clade, country, date, seq=True, ins="2999-01-01T00:00:00": {"accession": acc, "date_inserted": ins, "clade": clade, "isolation_country": country, "collection_date": date, "file_path": "x" if seq else None,
                                                                         "_qc": {}, "_excluded": False}
    rows = [g("OK", "Ib", "Uganda", "2024-05-01"), g("NOCLADE", "", "Uganda", "2024-05-01"), g("IIB", "IIb", "Uganda", "2024-05-01"), g("FOREIGN", "Ib", "Peru", "2024-05-01"),
            g("OLD", "Ib", "Uganda", "2019-05-01"), g("NOSEQ", "Ib", "Uganda", "2024-05-01", seq=False), g("BEFORE", "IIb", "Peru", "2024-05-01", ins="2000-01-01T00:00:00")]
    monkeypatch.setattr(export, "_load", lambda wid: ({}, rows))
    monkeypatch.setattr(export, "nextstrain_row", lambda r: {"clade": r["clade"], "country": r["isolation_country"], "date": r["collection_date"], "length": 190000, "PPX_accession": r["accession"]})
    monkeypatch.setattr(private, "rows_for_export", lambda wid: [])
    n = datasets.newcomers(ds, since="2001-01-01T00:00:00")
    assert n["new_total"] == 6 and n["members"] == 1                                         # BEFORE came in earlier; only OK matches every rule
    why = {r["id"]: r["why"] for r in n["rows"]}
    assert why == {"NOCLADE": "clade", "IIB": "clade", "FOREIGN": "country", "OLD": "year", "NOSEQ": "sequence"}
    a = datasets.add_manual(ds["id"], ["NOCLADE", "OLD", "NOT-A-GENOME"])
    assert a["added"] == 2
    ds2 = datasets.get(ds["id"])
    assert set(ds2["manual_include"]) == {"NOCLADE", "OLD"}
    n2 = datasets.newcomers(ds2, since="2001-01-01T00:00:00")
    assert n2["members"] == 3 and {r["id"] for r in n2["rows"]} == {"IIB", "FOREIGN", "NOSEQ"}      # added by hand: members whatever the rules say
    import pytest
    with pytest.raises(ValueError):
        datasets.add_manual(ds["id"], ["NOT-A-GENOME"])


def test_adaptive_candidates_need_persistence_spread_and_are_not_ancestral(monkeypatch):
    from collections import Counter
    from genomecurator import lineages, mutgantt
    monkeypatch.setattr(lineages, "chain", lambda ob, lin: [{"name": lin}])
    months = [f"2025-{m:02d}" for m in range(1, 13)]
    pool = Counter({m: 10 for m in months})                                       # 10 sampled genomes per month, 120 in all
    meta = {}
    cs_spread, cs_anc, cs_def = set(), set(), set()
    for i in range(120):
        m = months[i // 10]
        meta[f"g{i}"] = {"lineage": f"L{i % 4}", "country": f"C{i % 5}", "outbreak": "o", "clade": "Ib", "date": m + "-15"}
        if m >= "2025-07" and i % 3 != 0:                                       # rare early, common late: a spread over four lineages and five countries
            cs_spread.add(f"g{i}")
        cs_anc.add(f"g{i}")                                                       # in every genome
        if i % 4 == 0:
            cs_def.add(f"g{i}")                                                   # only L0: a lineage's own mutation
    def tr(cs, key, effect="nonsynonymous", gene="OPG1", apo=""):
        ms = [meta[s]["date"][:7] for s in cs]
        return mutgantt.adaptive_traits(key, cs, ms, pool, set(months[-6:]), meta, {(2, "G"): ["L0"]}, {"effect": effect, "gene": gene, "profile_of_interest": apo}, len(meta))
    a = tr(cs_spread, (1, "T"))
    assert a["tier"] in ("strong", "candidate") and a["rising"] and a["persistent"] and a["lineages"] == 4 and a["countries"] == 5 and a["score"] >= 6
    assert tr(cs_anc, (3, "T"))["tier"] == "fixed" and tr(cs_anc, (3, "T"))["score"] == 0                    # carried by every genome: ancestral, not an adaptation
    assert tr(cs_def, (2, "G"))["tier"] == "defining"                                                         # L0's own defining mutation
    assert tr(cs_spread, (1, "T"), effect="synonymous")["score"] == a["score"] - 1                           # no protein change: one point less
    assert tr(cs_spread, (1, "T"), apo="apobec")["tier"] == "apobec"                                          # APOBEC3-type: kept apart (a host edit, not evidence of selection)
    old = {f"g{i}" for i in range(0, 20)}                                                                     # seen only in the first two months: gone, not kept
    assert tr(old, (4, "A"))["persistent"] is False and tr(old, (4, "A"))["tier"] == ""


def test_co_travel_links_mutations_that_share_their_genomes():
    from genomecurator import mutgantt
    row = lambda pos, alt, first, cs: {"pos": pos, "alt": alt, "first": first, "_cs": set(cs)}
    rows = [row(10, "A", "2025-02", range(0, 10)), row(20, "C", "2025-01", range(0, 10)), row(30, "G", "2025-03", range(0, 6)), row(40, "T", "2025-04", range(50, 60))]
    links = mutgantt.co_travel(rows)
    by = {(x["a"], x["b"]): x for x in links}
    assert by[("20C", "10A")]["j"] == 1.0 and by[("20C", "10A")]["shared"] == 10                 # same genomes: a = the earlier-observed one
    assert ("30G", "10A") not in by and (("10A", "30G") in by and by[("10A", "30G")]["j"] == 0.6)  # 6 of 10 genomes: Jaccard 0.6 is the cut
    assert not any("40T" in k for k in by)                                                       # shares nothing: not linked


def test_independent_origins_count_settled_events_on_the_tree(monkeypatch):
    from genomecurator import lineages, mutprofile, nextstrain, origins, outbreak
    mk = lambda name, muts=(), kids=None: {"name": name, "branch_attrs": {"mutations": {"nuc": list(muts)}}, **({"children": kids} if kids else {})}
    # 100T arose on three branches: two with two descendants each, one with a single tip; 200G once on a branch with four tips
    tree = mk("root", [], [mk("n1", ["A100T", "C200G"], [mk("t1"), mk("t2"), mk("t3"), mk("t4")]), mk("n2", ["A100T"], [mk("t5"), mk("t6")]), mk("n3", [], [mk("t7", ["A100T"]), mk("t8", ["G300A"])])])
    monkeypatch.setattr(outbreak, "_tree", lambda run: (tree, "", {"auspice_file": "x.json", "build": "b"}))
    monkeypatch.setattr(outbreak, "_tree_reference", lambda ns: ("REF1", "A" * 400, []))
    monkeypatch.setattr(mutprofile, "summary", lambda run: {"reference": "REF1"})
    monkeypatch.setattr(lineages, "_genomes", lambda run: [{"id": f"t{i}", "lineage": "L1" if i < 5 else "L2", "country": "X" if i % 2 else "Y", "date": "2025-01-01", "outbreak": "", "clade": ""} for i in range(1, 9)])
    class F:
        def stat(self):
            class S: st_mtime = 1.0
            return S()
    monkeypatch.setattr(nextstrain, "auspice_dir", lambda: type("P", (), {"__truediv__": lambda s, o: F()})())
    origins._CACHE.clear()
    o = origins.independent_origins({"id": 1})
    assert {k: o[(100, "T")][k] for k in ("events", "settled", "singles", "tips", "lineages")} == {"events": 3, "settled": 2, "singles": 1, "tips": 6, "lineages": {"L1": 4, "L2": 2}}
    assert o[(200, "G")]["settled"] == 1 and o[(200, "G")]["tips"] == 4 and o[(300, "A")]["settled"] == 0 and o[(300, "A")]["singles"] == 1
    monkeypatch.setattr(mutprofile, "summary", lambda run: {"reference": "OTHER"})
    origins._CACHE.clear()
    assert origins.independent_origins({"id": 1}) is None                                  # a tree on another reference: no counts (the lineage approximation stays)
    # the evidence uses them: persistent + recurrent + rising + several countries + protein change -> a strong candidate; one origin -> no recurrence points
    from collections import Counter
    from genomecurator import mutgantt
    monkeypatch.setattr(lineages, "chain", lambda ob, lin: [{"name": lin}])
    months = [f"2025-{m:02d}" for m in range(1, 13)]
    pool = Counter({m: 10 for m in months})
    meta = {f"g{i}": {"lineage": f"L{i % 4}", "country": f"C{i % 5}", "outbreak": "o", "clade": "Ib", "date": months[i // 10] + "-15"} for i in range(120)}
    cs = {f"g{i}" for i in range(120) if months[i // 10] >= "2025-07" and i % 3 != 0}
    ms = [meta[s]["date"][:7] for s in cs]
    row = {"effect": "nonsynonymous", "gene": "OPG1", "profile_of_interest": ""}
    rec = mutgantt.adaptive_traits((1, "T"), cs, ms, pool, set(months[-6:]), meta, {}, row, 120, {"events": 3, "settled": 3, "singles": 0, "tips": 30}, True)
    one = mutgantt.adaptive_traits((1, "T"), cs, ms, pool, set(months[-6:]), meta, {}, row, 120, {"events": 1, "settled": 1, "singles": 0, "tips": 30}, True)
    assert rec["tier"] == "strong" and rec["origins"]["settled"] == 3 and one["tier"] != "strong" and rec["score"] - one["score"] == 3
    assert any("독립적으로 3번" in w for w in rec["why"]) and any("한 번만" in w for w in one["why"])


def test_neutral_baseline_judges_recurrence_against_synonymous_mutations(monkeypatch):
    from genomecurator import mutprofile, origins
    # 60 neutral (synonymous) transitions: 54 never settled twice... 6 settled once; protein-changing: one settled 4 times, one settled once, one never
    refseq = "A" * 3000
    o = {}
    for i in range(60):
        o[(100 + i, "G")] = {"events": 1, "settled": 1 if i < 6 else 0, "singles": 0 if i < 6 else 1, "tips": 2, "lineages": {}, "countries": {}}
    o[(1000, "G")] = {"events": 4, "settled": 4, "singles": 0, "tips": 12, "lineages": {}, "countries": {}}
    o[(1001, "G")] = {"events": 1, "settled": 1, "singles": 0, "tips": 3, "lineages": {}, "countries": {}}
    o[(1002, "G")] = {"events": 1, "settled": 0, "singles": 1, "tips": 0, "lineages": {}, "countries": {}}
    monkeypatch.setattr(mutprofile, "annotate", lambda pos, ref, alt, cds, seq: {"effect": "synonymous" if pos < 1000 else "nonsynonymous"})
    origins.neutral_baseline(o, refseq, [])
    hi, one, none = o[(1000, "G")], o[(1001, "G")], o[(1002, "G")]
    assert hi["class"] == "transition" and hi["neutral_n"] == 60
    assert hi["p_neutral"] < 1e-3 and hi["q_neutral"] < 0.05                                 # no neutral mutation reached 4 origins: surprising
    assert 0.05 < one["p_neutral"] < 0.2 and one["q_neutral"] > 0.05 and none["p_neutral"] == 1.0   # one origin is common among neutral mutations too
    assert abs(origins._count_sf(2, 0.4, 0.4) - 0.0616) < 1e-3 and origins._count_sf(2, 0.4, 1.0) > origins._count_sf(2, 0.4, 0.4)   # overdispersion widens the tail
    few = {(100 + i, "G"): {"events": 1, "settled": 0, "singles": 1, "tips": 0, "lineages": {}, "countries": {}} for i in range(5)}
    few[(1000, "G")] = {"events": 3, "settled": 3, "singles": 0, "tips": 9, "lineages": {}, "countries": {}}
    origins.neutral_baseline(few, refseq, [])
    assert few[(1000, "G")]["p_neutral"] is None                                              # fewer than 20 neutral mutations of the class: no p value
    assert o[(100, "G")]["p_neutral"] is None                                                # neutral mutations themselves are not tested


def test_trend_holds_the_country_mix_fixed():
    import random
    from genomecurator import trend
    rnd = random.Random(3)
    ids, month_of, country_of, carriers_mix, carriers_real = [], {}, {}, set(), set()
    for i in range(400):
        c = "A" if i < 200 else "B"                                        # A sampled in months 0-5, B in months 6-11
        m = rnd.randrange(0, 6) if c == "A" else rnd.randrange(6, 12)
        s = f"g{i}"
        ids.append(s)
        month_of[s], country_of[s] = m, c
        if rnd.random() < (0.1 if c == "A" else 0.7):                          # carrying depends on the country only: a rise that is just the country mix
            carriers_mix.add(s)
        if rnd.random() < 0.05 + 0.07 * m:                                      # carrying rises with time in both countries: a real rise
            carriers_real.add(s)
    ctx = trend.context(ids, month_of, country_of)
    mix = trend.carrier_trend(carriers_mix, ctx)
    real = trend.carrier_trend(carriers_real, ctx)
    assert mix["lo"] < 0 < mix["hi"]                                          # the country mix explains it: the interval spans 0
    assert real["lo"] > 0 and real["slope"] > 0.1
    assert trend.carrier_trend({"g0"}, ctx) is None                           # too few carriers


def test_hyphy_export_cuts_codons_in_coding_direction_and_writes_a_pruned_tree():
    from genomecurator import hyphyexport
    # alignment columns 1..12; a + gene on columns 2-7 (two codons), a - gene on columns 6-11 (read from the other strand)
    rows = {"g1": "AATGAAATTTCC", "g2": "AATGAC-NTTCC"}
    inv = {p: p for p in range(1, 13)}                                          # run reference position == alignment column here
    plus = {"gene": "P", "strand": "+", "exons": [(2, 7)]}
    minus = {"gene": "M", "strand": "-", "exons": [(6, 11)]}
    out, why = hyphyexport.codon_alignment(rows, inv, plus)
    assert out == {"g1": "ATGAAA", "g2": "ATGAC-"} and why == ""                                 # + strand: columns 2-7 as they are; the deleted base stays "-"
    om, _ = hyphyexport.codon_alignment(rows, inv, minus)
    assert om == {"g1": "GAAATT", "g2": "GAAN-G"}                                                  # - strand: columns 11..6 read backwards and complemented, N stays N
    assert hyphyexport.codon_alignment(rows, {p: p for p in range(1, 11)}, minus)[0] is None       # position 11 has no column: gene skipped
    assert hyphyexport.codon_alignment(rows, inv, {"gene": "B", "strand": "+", "exons": [(2, 6)]})[0] is None   # not a multiple of three
    assert hyphyexport.gene_positions(minus)[0] == (11, "-") and hyphyexport.gene_positions(minus)[-1] == (6, "-")
    mk = lambda name, div=None, kids=None: {"name": name, "node_attrs": {"div": div} if div is not None else {}, **({"children": kids} if kids else {})}
    tree = mk("root", 0, [mk("n1", 10, [mk("a", 12), mk("b", 15)]), mk("c", 4), mk("private:x", 3)])
    nw = hyphyexport.newick(tree, {"a", "b", "c"}, 100)
    assert nw.count("(") == 2 and "private:x" not in nw and nw.startswith("((a:0.02000000,b:0.05000000):0.10000000,c:0.04000000):0.00000000;")


def test_four_model_comparison_prepares_four_analyses_and_ranks_them(tmp_path, monkeypatch):
    import random
    import pytest
    from genomecurator import beastprep
    s = beastprep.parse_bets_settings({"compare": True})
    assert s["compare"] and [m["tag"] for m in beastprep._models(s)] == ["strict-constant", "strict-exponential", "relaxed-constant", "relaxed-exponential"] and all(m["dated"] for m in beastprep._models(s))
    assert [m["tag"] for m in beastprep._models(beastprep.parse_bets_settings({}))] == ["dated", "undated"]            # BETS keeps its two analyses
    rng = random.Random(2)
    names = [f"s{i}" for i in range(12)]
    seqs = {n: "".join(rng.choice("ACGT") for _ in range(60)) for n in names}
    real = {n: 2024 + i / 10 for i, n in enumerate(names)}
    monkeypatch.setattr(beastprep, "_bets_data", lambda cr, st: (seqs, real, 1e-5, 12, 1.1))
    cr = {"run_dir": str(tmp_path)}
    info = beastprep.bets_prepare(cr, {"compare": True, "estimator": "both"})
    root = tmp_path / "beast" / "bets" / "models"
    assert set(info["folders"]) == {"strict-constant", "strict-exponential", "relaxed-constant", "relaxed-exponential"} and (root / "settings.json").exists()
    x = {t: (root / t / "psss" / "beast.xml").read_text(encoding="utf-8") for t in info["folders"]}
    assert "ucld" in x["relaxed-constant"] and "ucld" not in x["strict-constant"]                                    # the clock differs
    assert x["strict-constant"] != x["strict-exponential"] and (root / "relaxed-exponential" / "mcmc" / "beast.xml").exists()   # so does the tree prior
    assert all(len({e.get("value") for e in __import__("xml.etree.ElementTree", fromlist=["x"]).fromstring(v).find("taxa").iter("date")}) > 1 for v in x.values())   # all with the real dates
    res = {"strict-constant": {"ps": -105.0, "ss": -105.4, "gss": -105.2}, "strict-exponential": {"ps": -101.0, "ss": -101.1, "gss": -101.0},
           "relaxed-constant": {"ps": -99.0, "ss": -99.2, "gss": -99.1}, "relaxed-exponential": {"ps": -90.0, "ss": -95.0, "gss": -93.0}}
    cmp = beastprep.compare_models(res)
    assert cmp["estimator"] == "gss" and [r["model"] for r in cmp["ranking"]] == ["relaxed-exponential", "relaxed-constant", "strict-exponential", "strict-constant"]
    assert cmp["ranking"][0]["log_bf_vs_best"] == 0 and cmp["ranking"][1]["log_bf_vs_best"] == -6.1 and cmp["ranking"][0]["estimator_spread"] == 5.0
    assert "relaxed-exponential" in cmp["text"] and "추정기들의 값이 2 이상" in cmp["text"]                           # the estimators disagree for the top model: said so
    assert beastprep.compare_models({"a": {"ps": None, "ss": None, "gss": None}})["ranking"] == []
