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
