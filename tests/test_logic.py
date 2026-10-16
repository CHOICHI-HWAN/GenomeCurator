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
