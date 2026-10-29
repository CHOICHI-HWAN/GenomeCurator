"""Cluster: queue rules, checksums, requeue on a silent worker, privacy routing, and a real worker against a real server (with a stand-in for BEAST)."""
import json
import socket
import stat
import threading
import time

import pytest


def _env(tmp_path, monkeypatch):
    from genomecurator import config
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    from genomecurator import cluster
    return cluster


def _pkg(tmp_path, name="beast1", tips=("A", "B")):
    d = tmp_path / "beast" / name
    d.mkdir(parents=True)
    (d / "beast.xml").write_text("<beast>toy</beast>", encoding="utf-8")
    (d / "settings.json").write_text(json.dumps({"seed": 7, "threads": 1, "chain_length": 1000, "tag": name, "tip_names": list(tips)}), encoding="utf-8")
    return d


def test_node_keys_are_hashed_and_revocable(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    n = c.create_node("w1")
    assert c.authenticate(n["key"])["name"] == "w1" and c.authenticate("wrong") is None and c.authenticate(None) is None
    from genomecurator.db import connect
    assert n["key"] not in str(dict(connect().execute("SELECT * FROM cluster_node").fetchone()))
    with pytest.raises(ValueError):
        c.create_node("w1")
    assert c.revoke_node(n["id"]) and c.authenticate(n["key"]) is None


def test_claim_order_one_job_per_worker_and_private_routing(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    con = connect()
    con.execute("INSERT INTO private_genome(watch_id,id) VALUES(1,'SECRET1')")
    con.commit()
    j_priv = c.submit_beast(_pkg(tmp_path, "p", ("A", "SECRET1")))
    j_pub = c.submit_beast(_pkg(tmp_path, "q"))
    assert c.get_job(j_priv)["private"] == 1 and c.get_job(j_pub)["private"] == 0
    plain = c.authenticate(c.create_node("plain")["key"])
    trusted = c.authenticate(c.create_node("trusted", allow_private=True)["key"])
    got = c.claim(plain)
    assert got["id"] == j_pub                       # the private job is skipped for a worker that does not accept it
    assert c.claim(plain) is None                   # one at a time
    assert c.claim(trusted)["id"] == j_priv
    assert got["inputs"][0]["name"] == "beast.xml" and len(got["inputs"][0]["sha256"]) == 64


def test_silent_worker_loses_the_job_and_attempts_run_out(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    jid = c.submit_beast(_pkg(tmp_path), max_attempts=2)
    n = c.authenticate(c.create_node("w")["key"])
    for attempt in (1, 2):
        assert c.claim(n)["id"] == jid
        con = connect()
        con.execute("UPDATE cluster_job SET heartbeat='2000-01-01T00:00:00', claimed='2000-01-01T00:00:00' WHERE id=?", (jid,))
        con.commit()
        assert c.reap() == 1
        assert c.get_job(jid)["state"] == ("queued" if attempt == 1 else "error")
    assert c.heartbeat(n, jid)["continue"] is False


def test_results_need_the_declared_checksum_and_all_required_files(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from genomecurator import cluster
    d = _pkg(tmp_path)
    jid = c.submit_beast(d)
    n = c.authenticate(c.create_node("w")["key"])
    other = c.authenticate(c.create_node("w2")["key"])
    c.claim(n)
    tmp = tmp_path / "up"
    tmp.write_text("log")
    with pytest.raises(PermissionError):
        c.store_output(other, jid, "beast.log", tmp, cluster.sha256(tmp))
    tmp.write_text("log")
    with pytest.raises(ValueError):
        c.store_output(n, jid, "../evil", tmp, cluster.sha256(tmp))      # names are plain file names only
    tmp.write_text("log")
    with pytest.raises(ValueError):
        c.store_output(n, jid, "beast.log", tmp, "0" * 64)
    assert not tmp.exists()                          # a rejected upload is deleted
    for name, text in (("beast.log", "log"),):
        tmp.write_text(text)
        c.store_output(n, jid, name, tmp, cluster.sha256(tmp))
    out = [{"name": "beast.log", "sha256": cluster.sha256(c.out_dir(jid) / "beast.log")}]
    r = c.complete(n, jid, out)
    assert r["ack"] is False and r["missing"] == ["beast.trees"]
    tmp.write_text("trees")
    c.store_output(n, jid, "beast.trees", tmp, cluster.sha256(tmp))
    out.append({"name": "beast.trees", "sha256": cluster.sha256(c.out_dir(jid) / "beast.trees")})
    assert c.complete(n, jid, out)["ack"] is True
    assert (d / "beast.log").read_text() == "log" and (d / "beast.trees").read_text() == "trees"
    assert json.loads((d / "cluster_job.json").read_text())["worker"] == "w"
    assert c.get_job(jid)["state"] == "done" and not (tmp_path / "cluster" / "jobs" / jid).exists()


def test_folder_must_be_inside_the_data_folder_and_cancel(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    out = tmp_path.parent / (tmp_path.name + "_outside")
    out.mkdir()
    (out / "beast.xml").write_text("x")
    with pytest.raises(ValueError):
        c.submit_beast(out)
    jid = c.submit_beast(_pkg(tmp_path))
    n = c.authenticate(c.create_node("w")["key"])
    c.claim(n)
    assert c.cancel(jid) and c.heartbeat(n, jid)["continue"] is False and c.claim(n) is None


def _serve(app):
    import uvicorn
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    srv = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=port, log_level="warning"))
    t = threading.Thread(target=srv.run, daemon=True)
    t.start()
    for _ in range(100):
        if srv.started:
            break
        time.sleep(0.05)
    return srv, port


def _fake_beast(tmp_path, fail=False):
    f = tmp_path / "fakebeast"
    f.write_text("#!/bin/sh\nif [ \"$1\" = \"-version\" ]; then echo 'BEAST v10.5.0'; exit 0; fi\n" + ("echo boom >&2; exit 3\n" if fail else
                 "printf 'state\\tposterior\\n0\\t-1\\n500\\t-1\\n1000\\t-1\\n' > beast.log\necho tree > beast.trees\necho ops > beast.ops\nsleep 1\n"))
    f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return str(f)


def test_a_worker_runs_a_job_end_to_end_and_cleans_up(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from genomecurator import web
    from genomecurator.clusterworker import Worker
    srv, port = _serve(web.create_app(scheduler=False))
    try:
        node = c.create_node("w1")
        d = _pkg(tmp_path)
        jid = c.submit_beast(d)
        w = Worker(f"http://127.0.0.1:{port}", node["key"], tmp_path / "wk", _fake_beast(tmp_path), poll=0.2, beat=0.5)
        t = threading.Thread(target=w.serve, daemon=True)
        t.start()
        for _ in range(200):
            if c.get_job(jid)["state"] in ("done", "error"):
                break
            time.sleep(0.1)
        w.stop.set()
        t.join(10)
        j = c.get_job(jid)
        assert j["state"] == "done", j
        assert (d / "beast.log").read_text().startswith("state") and (d / "beast.trees").read_text() == "tree\n"
        assert not any((tmp_path / "wk" / "jobs").glob("*"))              # the worker kept nothing
        assert any(n["name"] == "w1" and n["caps"]["tools"]["beast"] == "10.5.0" for n in c.list_nodes())
    finally:
        srv.should_exit = True


def test_a_failing_run_is_reported_and_retried_then_stops(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from genomecurator import web
    from genomecurator.clusterworker import Worker
    srv, port = _serve(web.create_app(scheduler=False))
    try:
        node = c.create_node("w1")
        jid = c.submit_beast(_pkg(tmp_path), max_attempts=2)
        w = Worker(f"http://127.0.0.1:{port}", node["key"], tmp_path / "wk", _fake_beast(tmp_path, fail=True), poll=0.2, beat=0.5)
        t = threading.Thread(target=w.serve, daemon=True)
        t.start()
        for _ in range(200):
            if c.get_job(jid)["state"] == "error":
                break
            time.sleep(0.1)
        w.stop.set()
        t.join(10)
        j = c.get_job(jid)
        assert j["state"] == "error" and j["attempts"] == 2 and "boom" in " ".join(x["msg"] for x in j["log"])
    finally:
        srv.should_exit = True


def test_results_wait_on_the_worker_while_the_main_server_is_away(tmp_path, monkeypatch):
    from genomecurator.clusterworker import Worker
    jd = tmp_path / "wk" / "jobs" / "abc"
    (jd / "run").mkdir(parents=True)
    (jd / "run" / "beast.log").write_text("x")
    (jd / "state.json").write_text(json.dumps({"phase": "finished", "job": {"id": "abc"}, "outputs": [{"name": "beast.log", "sha256": "0" * 64, "bytes": 1}]}))
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()                                         # nothing listens here
    w = Worker(f"http://127.0.0.1:{port}", "k", tmp_path / "wk", None)
    assert w.deliver(jd) is False and (jd / "run" / "beast.log").exists()


def test_http_api_needs_a_key_and_admin_is_local(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    cl = TestClient(web.create_app(scheduler=False))
    assert cl.post("/api/cluster/w/claim", json={}).status_code == 403
    assert cl.post("/api/cluster/w/claim", json={}, headers={"Authorization": "Bearer nope"}).status_code == 403
    r = cl.post("/api/cluster/nodes", json={"name": "n1"})
    assert r.status_code == 200 and "key" in r.json()
    assert cl.post("/api/cluster/nodes", json={"name": "n1"}).status_code == 400
    assert cl.get("/cluster").status_code == 200 and cl.get("/api/cluster/state").json()["nodes"][0]["name"] == "n1"
    assert cl.post("/api/cluster/w/claim", json={}, headers={"Authorization": "Bearer " + r.json()["key"]}).json() == {"job": None}


def test_other_machines_reach_the_worker_api_only(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    remote = TestClient(web.create_app(scheduler=False), client=("10.0.0.7", 5000))
    assert remote.get("/").status_code == 403 and remote.get("/api/cluster/state").status_code == 403
    assert remote.post("/api/cluster/nodes", json={"name": "x"}).status_code == 403
    assert remote.post("/api/cluster/w/claim", json={}).status_code == 403          # reaches the API, which wants a key
    assert "worker API" not in remote.post("/api/cluster/w/claim", json={}).text
    monkeypatch.setenv("GC_ALLOW_REMOTE_UI", "1")
    assert remote.get("/cluster").status_code == 200


def test_activating_a_server_as_a_sub_server_runs_jobs_and_survives_a_restart(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from genomecurator import clusterworker, web
    monkeypatch.setenv("PATH", str(tmp_path) + ":" + __import__("os").environ["PATH"])
    _fake_beast(tmp_path)                                   # found as "beast" on the PATH of this "sub-server"
    (tmp_path / "fakebeast").rename(tmp_path / "beastbin")
    monkeypatch.setattr(clusterworker, "find_beast", lambda explicit=None: str(tmp_path / "beastbin"))
    srv, port = _serve(web.create_app(scheduler=False))
    try:
        key = c.create_node("sub1")["key"]
        jid = c.submit_beast(_pkg(tmp_path))
        with pytest.raises(ValueError):
            clusterworker.sub_start("not-an-address", key)
        st = clusterworker.sub_start(f"http://127.0.0.1:{port}", key)
        assert st["enabled"] and st["running"]
        for _ in range(200):
            if c.get_job(jid)["state"] == "done":
                break
            time.sleep(0.1)
        assert c.get_job(jid)["state"] == "done"
        cfgf = tmp_path / "cluster" / "worker.json"
        assert (cfgf.stat().st_mode & 0o077) == 0                       # the key is readable by the owner only
        clusterworker.sub_stop(wait=True)
        assert not clusterworker.sub_status()["running"] and json.loads(cfgf.read_text())["enabled"] is False
        clusterworker.sub_start()                                        # the stored address and key are enough
        assert clusterworker.sub_status()["running"]
        clusterworker.sub_autostart()                                    # (re)starts only when enabled; harmless otherwise
        clusterworker.sub_stop(wait=True, forget=True)
        assert "key" not in json.loads(cfgf.read_text())
    finally:
        clusterworker.sub_stop(wait=True)
        srv.should_exit = True
