"""Cluster: queue rules, checksums, requeue on a silent worker, privacy routing, and a real worker against a real server (with a stand-in for BEAST)."""
import json
import socket
import stat
import threading
import time

import pytest

from genomecurator import clusterworker as _cw

_REAL_ENV_CHECK = _cw.env_check
OK_ENV = {"ok": True, "checks": [{"name": "BEAST", "ok": True, "detail": "BEAST 10.5.0"}, {"name": "BEAGLE", "ok": True, "detail": "library found"},
                                 {"name": "disk", "ok": True, "detail": "100 GB free"}], "beast_version": "10.5.0"}


@pytest.fixture(autouse=True)
def _environment_passes(tmp_path, monkeypatch):
    """The stand-in BEAST has no BEAGLE library next to it; the environment check itself is tested on its own below.
    Every test works in its own data folder: the clean-up below writes the worker settings, and must never reach the real ones."""
    from genomecurator import config
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DB_PATH", tmp_path / "t.db")
    monkeypatch.setattr(_cw, "env_check", lambda beast, workdir: OK_ENV)
    monkeypatch.setattr(_cw, "IN_PROCESS", True)           # the worker as a thread of the test (the separate process is tested on its own)
    _cw._env_cache.update({"at": 0.0, "env": None})
    yield
    _cw.sub_stop(wait=True)


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
    import datetime as dt
    c = _env(tmp_path, monkeypatch)
    monkeypatch.setattr(c, "_UP_SINCE", dt.datetime.now() - dt.timedelta(seconds=c.LEASE * 2))   # not a fresh start
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


def _fake_beast(tmp_path, fail=False, seconds=1):
    f = tmp_path / "fakebeast"
    f.write_text("#!/bin/sh\nif [ \"$1\" = \"-version\" ]; then echo 'BEAST v10.5.0'; exit 0; fi\n" + ("echo boom >&2; exit 3\n" if fail else
                 "printf 'state\\tposterior\\n0\\t-1\\n500\\t-1\\n' > beast.log\nsleep " + str(seconds) + "\nprintf '1000\\t-1\\n' >> beast.log\necho tree > beast.trees\necho ops > beast.ops\n"))
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
    remote = TestClient(web.create_app(scheduler=False), client=("8.8.8.8", 5000))
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


def test_same_network_machines_get_the_sub_server_page_only(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    app = web.create_app(scheduler=False)
    lan = TestClient(app, client=("192.168.0.20", 5000), follow_redirects=False)
    r = lan.get("/")
    assert r.status_code == 200 and "sub-server page" in r.text      # named, not silently redirected to this machine's page
    assert lan.get("/cluster/sub").status_code == 200
    s = lan.get("/api/cluster/sub").json()
    assert {"enabled", "running", "metrics", "install"} <= set(s) and s["metrics"]["cpus"] and "disk_free" in s["metrics"]
    for path in ("/genomes", "/private", "/api/cluster/state", "/cluster"):
        assert lan.get(path).status_code == 403
    assert lan.post("/api/cluster/nodes", json={"name": "x"}).status_code == 403              # keys are made on the main server
    assert lan.post("/api/cluster/sub/activate", json={"main": "bad", "key": "k"}).status_code == 400    # reaches the handler, which checks the address
    public = TestClient(app, client=("8.8.8.8", 5000), follow_redirects=False)
    assert public.get("/").status_code == 403 and public.get("/cluster/sub").status_code == 403
    monkeypatch.setenv("GC_SUB_PAGE_LAN", "0")
    assert lan.get("/cluster/sub").status_code == 403


def test_install_needs_a_package_manager_and_refuses_when_beast_exists(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    from genomecurator import beastprep, clusterworker
    monkeypatch.setattr(clusterworker.shutil, "which", lambda x: None)
    with pytest.raises(RuntimeError):
        clusterworker.install_beast()
    monkeypatch.setattr(clusterworker.shutil, "which", lambda x: "/bin/true")
    monkeypatch.setattr(beastprep, "status", lambda: {"ok": True, "version": "10.5.0", "how": ""})
    with pytest.raises(ValueError):
        clusterworker.install_beast()


def test_metrics_reach_the_main_server(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    n = c.authenticate(c.create_node("w")["key"])
    c.hello(n, {"kinds": ["beast"], "metrics": {"cpu_pct": 12.5, "ts": "2026-01-01T00:00:00"}})
    assert c.list_nodes()[0]["metrics"]["cpu_pct"] == 12.5
    c.claim(n, ["beast"], {"cpu_pct": 55.0, "ts": "x"})
    assert c.list_nodes()[0]["metrics"]["cpu_pct"] == 55.0



def test_environment_check_names_what_is_missing(tmp_path, monkeypatch):
    from genomecurator import beastprep
    monkeypatch.setattr(beastprep, "env_prefix", lambda: None)                      # do not find the BEAST of this machine
    monkeypatch.setattr(_cw, "MIN_FREE_GB", 0.0)
    e = _REAL_ENV_CHECK(None, tmp_path)
    assert not e["ok"] and {c["name"]: c["ok"] for c in e["checks"]} == {"BEAST": False, "BEAGLE": False, "disk": True}
    import shutil
    (tmp_path / "env" / "bin").mkdir(parents=True)
    fake = shutil.copy(_fake_beast(tmp_path), tmp_path / "env" / "bin" / "beast")       # a conda-like layout: <env>/bin/beast, <env>/lib
    (tmp_path / "env" / "lib").mkdir()
    e = _REAL_ENV_CHECK(fake, tmp_path)
    assert e["checks"][0]["ok"] and e["beast_version"] == "10.5.0" and not e["checks"][1]["ok"]      # BEAST starts, BEAGLE missing
    (tmp_path / "env" / "lib" / "libhmsbeagle-cpu.so").write_text("")
    assert _REAL_ENV_CHECK(fake, tmp_path)["ok"]
    monkeypatch.setattr(_cw, "MIN_FREE_GB", 1e9)
    assert not _REAL_ENV_CHECK(fake, tmp_path)["ok"]


def test_join_request_needs_approval_and_the_key_is_handed_over_once(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    with pytest.raises(ValueError):
        c.request_join("s1", "10.0.0.5", "x" * 30, {"env": {"ok": False}})          # a server whose environment failed may not ask
    rid = c.request_join("s1", "10.0.0.5", "s" * 30, {"env": OK_ENV, "cpus": 8})
    assert [r["name"] for r in c.list_requests()] == ["s1"] and "secret_sha" not in c.list_requests()[0]
    with pytest.raises(PermissionError):
        c.collect_request(rid, "wrong" * 6)
    assert c.collect_request(rid, "s" * 30) == {"state": "pending"} and c.list_nodes() == []
    assert c.decide_request(rid, True, allow_private=True)
    got = c.collect_request(rid, "s" * 30)
    assert got["state"] == "approved" and c.authenticate(got["key"])["allow_private"] == 1
    assert "key" not in c.collect_request(rid, "s" * 30)                            # only once
    rid2 = c.request_join("s2", "10.0.0.6", "t" * 30, {"env": OK_ENV})
    assert c.decide_request(rid2, False) and c.collect_request(rid2, "t" * 30)["state"] == "rejected"
    assert not c.decide_request(rid2, True)                                          # decided already


def test_a_sub_server_asks_waits_for_approval_and_starts_working(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    monkeypatch.setattr(_cw, "find_beast", lambda explicit=None: _fake_beast(tmp_path))
    srv, port = _serve(web.create_app(scheduler=False))
    try:
        jid = c.submit_beast(_pkg(tmp_path))
        st = _cw.sub_request(f"http://127.0.0.1:{port}", name="sub-a")
        assert st["join"]["state"] == "pending" and not st["running"]
        time.sleep(0.5)
        assert c.get_job(jid)["state"] == "queued"                                  # nothing runs before approval
        main = TestClient(web.create_app(scheduler=False))
        req = main.get("/api/cluster/state").json()["requests"]
        assert req[0]["name"] == "sub-a" and req[0]["info"]["env"]["ok"]
        assert main.post(f"/api/cluster/requests/{req[0]['id']}/approve", json={}).json()["ok"]
        for _ in range(300):
            if c.get_job(jid)["state"] == "done":
                break
            time.sleep(0.1)
        assert c.get_job(jid)["state"] == "done"
        node = c.list_nodes()[0]
        assert node["name"] == "sub-a" and node["caps"]["env"]["ok"] and node["caps"]["app_version"]
        s = _cw.sub_status()
        assert s["running"] and s["join"]["state"] == "approved" and s["history"][0]["result"] == "done"
        assert main.post(f"/api/cluster/nodes/{node['id']}/recheck").json()["ok"]
    finally:
        _cw.sub_stop(wait=True)
        srv.should_exit = True


def test_a_sub_server_with_a_broken_environment_cannot_ask_or_start(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    bad = {"ok": False, "checks": [{"name": "BEAGLE", "ok": False, "detail": "missing"}]}
    monkeypatch.setattr(_cw, "env_check", lambda beast, workdir: bad)
    with pytest.raises(ValueError, match="BEAGLE"):
        _cw.sub_request("http://127.0.0.1:9")
    with pytest.raises(ValueError, match="BEAGLE"):
        _cw.sub_start("http://127.0.0.1:9", "k" * 30)



def test_a_silent_job_shows_its_grace_period_before_it_is_handed_out_again(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from datetime import datetime, timedelta
    from genomecurator.db import connect
    jid = c.submit_beast(_pkg(tmp_path))
    n = c.authenticate(c.create_node("w")["key"])
    c.claim(n)
    assert c.get_job(jid)["stale"] is False
    t = (datetime.now() - timedelta(seconds=c.STALE + 5)).isoformat(timespec="seconds")
    con = connect()
    con.execute("UPDATE cluster_job SET heartbeat=? WHERE id=?", (t, jid))
    con.commit()
    j = c.get_job(jid)
    assert j["stale"] and 0 < j["requeue_in"] <= c.LEASE - c.STALE and c.reap() == 0        # in its grace period: still running
    assert c.list_nodes()[0]["running"]["stale"]


def test_stop_on_the_sub_server_gives_the_job_back_at_once_without_using_an_attempt(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from genomecurator import web
    monkeypatch.setattr(_cw, "find_beast", lambda explicit=None: _fake_beast(tmp_path, seconds=30))
    srv, port = _serve(web.create_app(scheduler=False))
    try:
        key = c.create_node("s")["key"]
        jid = c.submit_beast(_pkg(tmp_path))
        _cw.sub_start(f"http://127.0.0.1:{port}", key)
        for _ in range(100):
            if c.get_job(jid)["state"] == "running":
                break
            time.sleep(0.1)
        assert c.get_job(jid)["state"] == "running"
        _cw.sub_pause(True)
        for _ in range(100):
            if c.get_job(jid)["state"] == "queued":
                break
            time.sleep(0.1)
        j = c.get_job(jid)
        assert j["state"] == "queued" and j["attempts"] == 0 and "not counted" in j["log"][-1]["msg"]
        time.sleep(1)
        assert c.get_job(jid)["state"] == "queued"                                # a stopped sub-server takes no new work
        assert _cw.sub_status()["paused"] and c.list_nodes()[0]["paused"]
        assert _cw.load_config()["paused"] is True                                 # kept across restarts
        _cw.sub_pause(False)
        for _ in range(100):
            if c.get_job(jid)["state"] == "running":
                break
            time.sleep(0.1)
        assert c.get_job(jid)["state"] == "running"
    finally:
        _cw.sub_stop(wait=True)
        srv.should_exit = True


def test_terminal_tokens_are_signed_short_lived_and_single_use(tmp_path, monkeypatch):
    import hashlib
    from genomecurator import clustershell as sh
    key = "k" * 43
    tok = sh.make_token(hashlib.sha256(key.encode()).hexdigest(), "n1")
    assert sh.check_token(tok, key)["n"] == "n1"
    with pytest.raises(ValueError, match="used"):
        sh.check_token(tok, key)
    with pytest.raises(ValueError, match="not valid"):
        sh.check_token(sh.make_token(hashlib.sha256(b"other").hexdigest(), "n1"), key)
    with pytest.raises(ValueError):
        sh.check_token(tok, None)
    monkeypatch.setattr(sh, "TOKEN_SECONDS", -1)
    with pytest.raises(ValueError, match="expired"):
        sh.check_token(sh.make_token(hashlib.sha256(key.encode()).hexdigest(), "n1"), key)


@pytest.mark.skipif(__import__("os").name != "posix", reason="pseudo-terminals")
def test_the_main_server_opens_a_terminal_on_a_sub_server_that_allows_it(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    node = c.create_node("sub1")
    n = c.authenticate(node["key"])
    c.hello(n, {"kinds": ["beast"], "ip": "127.0.0.1", "port": "8011", "metrics": {"shell": False, "ts": "x"}})
    main = TestClient(web.create_app(scheduler=False))
    assert main.post(f"/api/cluster/nodes/{node['id']}/shell").status_code == 409            # not allowed by the sub-server
    c.claim(n, [], {"shell": True, "ts": "y"})
    r = main.post(f"/api/cluster/nodes/{node['id']}/shell").json()
    assert r["url"] == "ws://127.0.0.1:8011/api/cluster/shell" and r["token"]
    assert TestClient(web.create_app(scheduler=False), client=("192.168.0.9", 1)).post(f"/api/cluster/nodes/{node['id']}/shell").status_code == 403
    sub = TestClient(web.create_app(scheduler=False))                                    # the same program on the sub-server
    _cw.save_config({"main": "http://x", "key": node["key"], "shell": False})
    with sub.websocket_connect("/api/cluster/shell") as ws:
        ws.send_text(json.dumps({"token": r["token"]}))
        assert "switched off" in ws.receive_text()
    _cw.sub_shell(True)
    tok = main.post(f"/api/cluster/nodes/{node['id']}/shell").json()["token"]
    with sub.websocket_connect("/api/cluster/shell") as ws:
        ws.send_text(json.dumps({"token": tok, "cols": 80, "rows": 24}))
        ws.send_text(json.dumps({"t": "i", "d": "echo gc-$((40+2))\n"}))
        seen = ""
        for _ in range(200):
            seen += ws.receive_text()
            if "gc-42" in seen:
                break
        assert "gc-42" in seen
        ws.send_text(json.dumps({"t": "i", "d": "exit\n"}))
    with sub.websocket_connect("/api/cluster/shell") as ws:                              # the same token cannot be used twice
        ws.send_text(json.dumps({"token": tok}))
        assert "used already" in ws.receive_text()


def _fake_netdata(new_names=True, claimed=False):
    """A stand-in Netdata agent: /api/v1/info and /api/v1/data for the four charts (also under /host/<name>/ like a streaming parent)."""
    import http.server
    now = int(time.time())
    charts = {"system.cpu": (["time", "user", "system"], lambda i: [60.0, 5.0 + i]),
              "system.ram": (["time", "free", "used", "cached", "buffers"], lambda i: [50.0, 30.0, 15.0, 5.0]),
              "system.load": (["time", "load1", "load5", "load15"], lambda i: [3.0, 2.0, 1.0]),
              "system.io": ((["time", "reads", "writes"] if new_names else ["time", "in", "out"]), lambda i: [2048.0, -4096.0]),
              "system.net": (["time", "received", "sent"], lambda i: [1000.0, -500.0])}

    class H(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def do_GET(self):
            from urllib.parse import urlparse, parse_qs
            u = urlparse(self.path)
            path = u.path.split("/api/", 1)[-1]
            if path == "v1/info":
                body = {"version": "v2.0.0", "agent-claimed": claimed}
            elif path == "v1/data":
                lab, fn = charts[parse_qs(u.query)["chart"][0]]
                body = {"labels": lab, "data": [[now - 60 * i] + fn(i) for i in range(5)]}           # newest first, like Netdata
            else:
                self.send_response(404)
                self.end_headers()
                return
            b = json.dumps(body).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(b)
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f"http://127.0.0.1:{srv.server_address[1]}"


@pytest.mark.parametrize("new_names", [True, False])
def test_netdata_history_is_read_oldest_first_in_plain_units(new_names):
    from genomecurator import netdata
    srv, base = _fake_netdata(new_names, claimed=True)
    try:
        assert netdata.probe(base) == {"ok": True, "version": "v2.0.0", "claimed": True, "port": srv.server_address[1]}
        h = netdata.history(base, minutes=10)
        assert h["t"] == sorted(h["t"]) and len(h["t"]) == 5
        assert h["cpu"][0] == 69.0 and h["cpu"][-1] == 65.0                     # oldest (i=4) first
        assert h["mem"][0] == 30.0 and h["load1"][0] == 3.0
        assert h["disk_read"][0] == 2.0 and h["disk_write"][0] == 4.0           # KiB/s -> MiB/s, sign dropped
        assert h["net_in"][0] == 1.0 and h["net_out"][0] == 0.5                # kbit/s -> Mbit/s
        assert netdata.history(base, minutes=10, host="sub-a")["source"].endswith("/host/sub-a")
    finally:
        srv.shutdown()
    assert netdata.probe("http://127.0.0.1:9") == {"ok": False}
    assert netdata.history_or_none("http://127.0.0.1:9")["ok"] is False


def test_main_server_shows_a_sub_servers_history_with_its_jobs(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import netdata, web
    srv, base = _fake_netdata()
    try:
        node = c.create_node("sub1")
        n = c.authenticate(node["key"])
        c.hello(n, {"kinds": ["beast"], "ip": "127.0.0.1", "host": "sub1", "netdata": {"ok": True, "port": srv.server_address[1]}})
        jid = c.submit_beast(_pkg(tmp_path))
        c.claim(n)
        main = TestClient(web.create_app(scheduler=False))
        h = main.get(f"/api/cluster/nodes/{node['id']}/history?minutes=30").json()
        assert h["ok"] and len(h["cpu"]) == 5 and h["jobs"][0]["label"].startswith("BEAST") and h["jobs"][0]["end"] is None
        assert main.get("/static/gc/minichart.js").status_code == 200
        # no agent on the sub-server: read it through the main server's own Netdata (a streaming parent), by host name
        c.hello(n, {"kinds": ["beast"], "ip": "127.0.0.1", "host": "sub1", "netdata": {"ok": False}})
        monkeypatch.setattr(netdata, "LOCAL", base)
        h = main.get(f"/api/cluster/nodes/{node['id']}/history").json()
        assert h["ok"] and h["source"].endswith("/host/sub1")
        sub = TestClient(web.create_app(scheduler=False))
        assert sub.get("/api/cluster/sub/history?minutes=30").json()["ok"]
        assert jid
    finally:
        srv.shutdown()


@pytest.mark.skipif(__import__("shutil").which("node") is None, reason="node is needed to check the page scripts")
def test_the_cluster_pages_scripts_parse(tmp_path, monkeypatch):
    """The pages are mostly JavaScript built from strings; a quote out of place breaks a whole page and no API test notices."""
    import re
    import subprocess
    c = _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    node = c.create_node("n1")
    cl = TestClient(web.create_app(scheduler=False))
    for path in ("/cluster", "/cluster/sub", f"/cluster/shell/{node['id']}"):
        html = cl.get(path).text
        js = "\n".join(re.findall(r"<script>(.*?)</script>", html, re.S))
        f = tmp_path / "page.js"
        f.write_text(js, encoding="utf-8")
        r = subprocess.run(["node", "--check", str(f)], capture_output=True, text=True)
        assert r.returncode == 0, f"{path}: {r.stderr[:400]}"
    from pathlib import Path
    r = subprocess.run(["node", "--check", str(Path(web.__file__).parent / "static" / "gc" / "minichart.js")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_netdata_streaming_guide_names_the_link_address_and_who_streams(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import netdata, web
    srv, base = _fake_netdata()
    try:
        monkeypatch.setattr(netdata, "LOCAL", base)
        monkeypatch.setattr(netdata, "streaming_hosts", lambda b=None: ["main", "sub1"])
        n = c.authenticate(c.create_node("sub1")["key"])
        c.hello(n, {"ip": "127.0.0.1", "host": "sub1"})
        j = TestClient(web.create_app(scheduler=False)).get("/api/cluster/netdata").json()
        key = netdata.stream_key()
        assert j["parent_address"] == "127.0.0.1" and f"[{key}]" in j["parent_conf"] and "allow from = 127.0.0.*" in j["parent_conf"]
        assert f"api key = {key}" in j["child_conf"] and "destination = 127.0.0.1:19999" in j["child_conf"]
        assert j["streaming"] == [{"name": "sub1", "host": "sub1", "streams": True}]
        assert netdata.stream_key() == key                                       # made once
    finally:
        srv.shutdown()


def test_a_browser_on_the_server_machine_is_not_treated_as_a_visitor(tmp_path, monkeypatch):
    """`serve --host 10.0.0.1` opened at http://10.0.0.1:8001 on that same machine: the client address equals the address it was sent to."""
    _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    app = web.create_app(scheduler=False)
    own = TestClient(app, client=("10.0.0.1", 5000), headers={"host": "10.0.0.1:8001"}, follow_redirects=False)
    assert own.get("/cluster").status_code == 200                                  # the main screens stay open on the server itself
    assert own.post("/api/cluster/nodes", json={"name": "n1"}).status_code == 200   # and so does administration
    other = TestClient(app, client=("10.0.0.2", 5000), headers={"host": "10.0.0.1:8001"}, follow_redirects=False)
    assert other.get("/cluster").status_code == 403
    r = other.get("/")
    assert r.status_code == 200 and "sub-server page" in r.text                     # named, not redirected to this machine's sub-server page
    assert "10.0.0.2" in r.text and "location" not in r.headers
    assert other.post("/api/cluster/nodes", json={"name": "n2"}).status_code == 403


def test_the_sub_server_page_says_which_machine_it_is(tmp_path, monkeypatch):
    import socket
    c = _env(tmp_path, monkeypatch)
    from fastapi.testclient import TestClient
    from genomecurator import web
    cl = TestClient(web.create_app(scheduler=False))
    j = cl.get("/api/cluster/sub").json()
    assert j["host"] == socket.gethostname() and j["workers_here"] == 0
    c.create_node("w1")
    assert cl.get("/api/cluster/sub").json()["workers_here"] == 1                   # the page can warn that this machine is also a main server


def test_the_main_server_works_as_a_worker_of_itself_by_default(tmp_path, monkeypatch):
    """It holds the data, so there is no reason for it to idle; the person can still revoke it or switch it off."""
    c = _env(tmp_path, monkeypatch)
    import socket
    from genomecurator.db import connect
    monkeypatch.setenv("GC_SERVE_PORT", "8001")
    monkeypatch.setattr(_cw, "sub_start", lambda *a, **k: {"started": True})      # no real worker thread in the test
    _cw.main_worker_autostart()
    nodes = c.list_nodes()
    assert [n["name"] for n in nodes] == [socket.gethostname()] and nodes[0]["allow_private"] == 1
    cfg = _cw.load_config()
    assert cfg["enabled"] and cfg["main"] == "http://127.0.0.1:8001" and c.authenticate(cfg["key"])["id"] == nodes[0]["id"]
    _cw.main_worker_autostart()
    assert len(c.list_nodes()) == 1                                               # not created again on the next start
    c.revoke_node(nodes[0]["id"])
    _cw.main_worker_autostart()
    assert c.list_nodes() == []                                                   # a revoked local worker stays revoked
    # an address of this machine also points here, not only 127.0.0.1
    assert _cw._points_here("http://127.0.0.1:8001", "8001") and _cw._points_here("http://localhost:8001", "8001")
    assert not _cw._points_here("http://10.0.0.251:8001", "8001") and not _cw._points_here("http://127.0.0.1:9999", "8001")
    # a machine that serves another main server is left alone
    _cw.save_config({"main": "http://10.0.0.251:8001", "key": "k" * 20, "enabled": True})
    con = connect(); con.execute("DELETE FROM app_meta WHERE key='cluster_local_worker'"); con.commit()
    _cw.main_worker_autostart()
    assert c.list_nodes() == [] and _cw.load_config()["main"] == "http://10.0.0.251:8001"
    monkeypatch.setenv("GC_MAIN_WORKER", "0")
    _cw.save_config({})
    _cw.main_worker_autostart()
    assert c.list_nodes() == []


def test_a_local_worker_whose_key_no_longer_works_is_registered_again(tmp_path, monkeypatch):
    """Leftover settings from an earlier cluster (switched off, key revoked) must not keep the main server idle for ever."""
    c = _env(tmp_path, monkeypatch)
    import socket
    monkeypatch.setenv("GC_SERVE_PORT", "8001")
    monkeypatch.setattr(_cw, "sub_start", lambda *a, **k: {"started": True})
    _cw.save_config({"main": "http://127.0.0.1:8001", "key": "dead" * 8, "enabled": False, "name": "old"})
    _cw.main_worker_autostart()
    assert [n["name"] for n in c.list_nodes()] == [socket.gethostname()]
    assert c.authenticate(_cw.load_config()["key"]) is not None and _cw.load_config()["enabled"] is True


def test_a_worker_fills_its_cores_with_several_jobs(tmp_path, monkeypatch):
    """One job per machine wastes a 32-core server; the main server hands out work until a worker's cores are taken."""
    c = _env(tmp_path, monkeypatch)
    n = c.authenticate(c.create_node("big")["key"])
    c.hello(n, {"kinds": ["beast"], "cpus": 32})
    ids = []
    for i in range(4):
        d = _pkg(tmp_path, f"j{i}")
        (d / "settings.json").write_text(json.dumps({"seed": 1, "threads": 8, "chain_length": 1000, "tag": f"j{i}", "tip_names": []}))
        ids.append(c.submit_beast(d))
    got = [c.claim(n) for _ in range(4)]
    assert [g["id"] for g in got if g] == ids[:3]            # 32 cores - 4 reserved = 28, three jobs of 8 threads fit
    assert got[3] is None and c.get_job(ids[3])["state"] == "queued"
    node = c.list_nodes()[0]
    assert len(node["running_jobs"]) == 3 and node["cpus_busy"] == 24 and node["cpus_free"] == 4
    small = c.authenticate(c.create_node("small")["key"])
    c.hello(small, {"kinds": ["beast"], "cpus": 8})
    assert c.claim(small)["id"] == ids[3]                     # the last one goes to the other worker
    assert c.claim(small) is None


def test_a_run_that_finished_while_the_worker_was_away_is_delivered_not_repeated(tmp_path, monkeypatch):
    """The run is its own process: restarting the worker must not throw away a chain that completed meanwhile."""
    c = _env(tmp_path, monkeypatch)
    from genomecurator import web
    srv, port = _serve(web.create_app(scheduler=False))
    try:
        key = c.create_node("w")["key"]
        d = _pkg(tmp_path)
        jid = c.submit_beast(d)
        n = c.authenticate(key)
        job = c.claim(n)                                       # the job is out with the worker
        w = _cw.Worker(f"http://127.0.0.1:{port}", key, tmp_path / "wk", None)
        jd = tmp_path / "wk" / "jobs" / jid
        (jd / "run").mkdir(parents=True)
        (jd / "run" / "beast.log").write_text("state\tposterior\n0\t-1\n1000\t-1\n")      # chain_length of _pkg is 1000
        for f in ("beast.trees", "beast.ops", "beast.screen.log"):
            (jd / "run" / f).write_text(f)
        (jd / "state.json").write_text(json.dumps({"phase": "running", "job": job}))
        w.resume()
        assert c.get_job(jid)["state"] == "done" and not jd.exists()
        assert (d / "beast.trees").read_text() == "beast.trees"
        # a chain that did not reach its length is not salvaged
        jid2 = c.submit_beast(_pkg(tmp_path, "half"))
        job2 = c.claim(n)
        jd2 = tmp_path / "wk" / "jobs" / jid2
        (jd2 / "run").mkdir(parents=True)
        (jd2 / "run" / "beast.log").write_text("state\tposterior\n0\t-1\n400\t-1\n")
        for f in ("beast.trees", "beast.ops", "beast.screen.log"):
            (jd2 / "run" / f).write_text(f)
        (jd2 / "state.json").write_text(json.dumps({"phase": "running", "job": job2}))
        w.resume()
        assert not jd2.exists() and c.get_job(jid2)["state"] == "running"     # left to the main server's grace period
    finally:
        srv.should_exit = True



def test_a_worker_started_from_the_app_is_a_separate_process_found_again_after_a_restart(tmp_path, monkeypatch):
    """Restarting the app must not stop the computing: the worker runs as its own process, once per data folder."""
    _env(tmp_path, monkeypatch)
    monkeypatch.setattr(_cw, "IN_PROCESS", False)
    started = []

    class P:
        pid = 424242
    monkeypatch.setattr(_cw.subprocess, "Popen", lambda cmd, **kw: started.append((cmd, kw)) or P())
    alive = {"pid": None}
    monkeypatch.setattr(_cw, "_worker_pid", lambda: alive["pid"])
    cfg = {"main": "http://127.0.0.1:8001", "key": "k" * 30, "enabled": True}
    _cw.save_config(cfg)
    _cw._start_worker(cfg)
    assert len(started) == 1 and started[0][0][-2:] == ["worker", "--from-config"]
    assert started[0][1].get("start_new_session") is True                     # its own session: it outlives the app
    assert _cw._pid_path().read_text() == "424242"
    alive["pid"] = 424242                                                      # the app restarts; the worker is still there
    _cw._start_worker(cfg)
    assert len(started) == 1                                                   # found again, not started twice
    _cw._status_path().write_text(json.dumps({"pid": 424242, "status": {"state": "running"}, "jobs": [{"id": "j1", "pct": 40.0}], "history": [], "beast": "b"}))
    st = _cw.sub_status()
    assert st["running"] and st["jobs"][0]["pct"] == 40.0 and st["status"]["state"] == "running"


def test_the_worker_process_follows_switch_off_and_pause_from_the_app(tmp_path, monkeypatch):
    _env(tmp_path, monkeypatch)
    w = _cw.Worker("http://127.0.0.1:9", "key1", tmp_path / "wk", None)
    _cw._apply_config(w, {"enabled": True, "key": "key1", "paused": True}, "key1")
    assert w.paused.is_set() and not w.stop.is_set()
    _cw._apply_config(w, {"enabled": True, "key": "key1", "paused": False}, "key1")
    assert not w.paused.is_set()
    _cw._apply_config(w, {"enabled": True, "key": "another", "paused": False}, "key1")
    assert w.stop.is_set()                                                     # a new key means a new worker
    w2 = _cw.Worker("http://127.0.0.1:9", "key1", tmp_path / "wk2", None)
    _cw._apply_config(w2, {"enabled": False, "key": "key1"}, "key1")
    assert w2.stop.is_set()
    _cw._write_status(w2)
    assert json.loads(_cw._status_path().read_text())["status"]["state"] == "starting"


def test_restarting_the_app_keeps_a_running_worker_process(tmp_path, monkeypatch):
    """The app's start-up calls sub_start(); with a worker process already running on the same settings it must not stop it."""
    _env(tmp_path, monkeypatch)
    monkeypatch.setattr(_cw, "IN_PROCESS", False)
    monkeypatch.setattr(_cw, "_worker_pid", lambda: 31337)                      # a worker process is alive
    stopped, spawned = [], []
    monkeypatch.setattr(_cw, "sub_stop", lambda **k: stopped.append(k) or {})
    monkeypatch.setattr(_cw, "_spawn_worker", lambda: spawned.append(1) or 1)
    _cw.save_config({"main": "http://127.0.0.1:8001", "key": "k" * 30, "enabled": True})
    _cw.sub_start()
    assert stopped == [] and spawned == []                                     # left alone
    _cw.sub_start("http://10.0.0.9:8001", "z" * 30)                            # other settings: then it is restarted
    assert stopped and _cw.load_config()["main"] == "http://10.0.0.9:8001"


def test_a_restarted_worker_hands_back_jobs_it_no_longer_runs(tmp_path, monkeypatch):
    c = _env(tmp_path, monkeypatch)
    n = c.authenticate(c.create_node("w")["key"])
    c.hello(n, {"kinds": ["beast"], "cpus": 64})
    a = c.submit_beast(_pkg(tmp_path, "a"))
    b = c.submit_beast(_pkg(tmp_path, "b"))
    c.claim(n); c.claim(n)
    r = c.sync_running(n, [a])                                                 # the worker still runs a, lost b
    assert r["released"] == [b]
    assert c.get_job(a)["state"] == "running" and c.get_job(b)["state"] == "queued" and c.get_job(b)["attempts"] == 0


def test_a_job_is_not_handed_out_again_just_because_the_main_server_restarted(tmp_path, monkeypatch):
    """The worker keeps computing while the main server is away; it must get its full lease to report in before the job is taken from it."""
    import datetime as dt
    c = _env(tmp_path, monkeypatch)
    from genomecurator.db import connect
    jid = c.submit_beast(_pkg(tmp_path))
    n = c.authenticate(c.create_node("w")["key"])
    c.claim(n)
    old = (dt.datetime.now() - dt.timedelta(seconds=c.LEASE * 3)).isoformat(timespec="seconds")
    con = connect()
    con.execute("UPDATE cluster_job SET heartbeat=?, claimed=? WHERE id=?", (old, old, jid))
    con.commit()
    monkeypatch.setattr(c, "_UP_SINCE", dt.datetime.now())                     # as if the server had just come back
    assert c.reap() == 0 and c.get_job(jid)["state"] == "running"
    monkeypatch.setattr(c, "_UP_SINCE", dt.datetime.now() - dt.timedelta(seconds=c.LEASE * 2))
    assert c.reap() == 1 and c.get_job(jid)["state"] == "queued"               # once it has been up long enough, the silence counts
