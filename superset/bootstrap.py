"""Create the GenomeCurator database connection, dataset, charts and dashboard in a running Superset
through its REST API (idempotent: existing objects are reused). Usage:
  python superset/bootstrap.py --url http://127.0.0.1:8089 --user admin --password-file superset/home/admin_password.txt --db data/genomecurator.db
"""
import argparse
import json
import sys
from pathlib import Path

import requests

VIEW = "v_genome_analysis"
DB_NAME = "GenomeCurator"
DASH = "GenomeCurator 메타데이터 분석"


def metric(label, col=None, agg="COUNT"):
    if col is None:
        return {"expressionType": "SQL", "sqlExpression": "COUNT(*)", "label": label}
    return {"expressionType": "SIMPLE", "column": {"column_name": col}, "aggregate": agg, "label": label}


COUNT = metric("유전체 수")
CHARTS = [
    # name, viz_type, params
    ("총 유전체 수", "big_number_total", {"metric": COUNT}),
    ("바이러스별 유전체 수", "pie", {"groupby": ["virus"], "metric": COUNT, "show_labels": True, "label_type": "key_value"}),
    ("수집 연도별 추이 (바이러스별)", "echarts_timeseries_bar",
     {"x_axis": "collection_year", "groupby": ["virus"], "metrics": [COUNT], "stack": True, "x_axis_sort_asc": True}),
    ("국가별 상위 20 (바이러스별)", "echarts_timeseries_bar",
     {"x_axis": "country", "groupby": ["virus"], "metrics": [COUNT], "stack": True, "orientation": "horizontal",
      "row_limit": 20, "x_axis_sort": "유전체 수", "x_axis_sort_asc": False}),
    ("숙주 × 국가 (히트맵)", "heatmap_v2",
     {"x_axis": "host", "groupby": "country", "metric": COUNT, "row_limit": 500, "normalize_across": "heatmap"}),
    ("유전체 길이 구간 × 상태", "echarts_timeseries_bar",
     {"x_axis": "length_bin", "groupby": ["genome_status"], "metrics": [COUNT], "stack": True}),
    ("바이러스 → 숙주 → 국가 (트리맵)", "treemap_v2", {"groupby": ["virus", "host", "country"], "metric": COUNT}),
    ("메타데이터 완성도: SRA / FASTA / Pathoplexus", "table",
     {"query_mode": "aggregate", "groupby": ["virus", "has_sra", "has_fasta", "in_pathoplexus"], "metrics": [COUNT], "row_limit": 1000}),
    ("Clade / Lineage 분포", "echarts_timeseries_bar",
     {"x_axis": "clade", "groupby": ["virus"], "metrics": [COUNT], "stack": True}),
    ("피벗: 수집 연도 × 국가", "pivot_table_v2",
     {"groupbyRows": ["country"], "groupbyColumns": ["collection_year"], "metrics": [COUNT], "aggregateFunction": "Sum",
      "rowTotals": True, "colTotals": True, "row_limit": 50000}),
]


class Api:
    def __init__(self, url, user, password):
        self.u = url.rstrip("/")
        self.s = requests.Session()
        r = self.s.post(f"{self.u}/api/v1/security/login", json={"username": user, "password": password, "provider": "db", "refresh": True})
        r.raise_for_status()
        self.s.headers["Authorization"] = "Bearer " + r.json()["access_token"]
        c = self.s.get(f"{self.u}/api/v1/security/csrf_token/")
        c.raise_for_status()
        self.s.headers.update({"X-CSRFToken": c.json()["result"], "Referer": self.u})

    def call(self, method, path, **kw):
        r = self.s.request(method, f"{self.u}/api/v1/{path}", **kw)
        if not r.ok:
            raise SystemExit(f"{method} {path} -> {r.status_code}: {r.text[:400]}")
        return r.json()

    def find(self, path, col, val):
        q = json.dumps({"filters": [{"col": col, "opr": "eq", "value": val}]})
        res = self.call("GET", f"{path}/", params={"q": q})["result"]
        return res[0]["id"] if res else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default="http://127.0.0.1:8089")
    ap.add_argument("--user", default="admin")
    ap.add_argument("--password-file", required=True)
    ap.add_argument("--db", required=True, help="GenomeCurator SQLite file")
    a = ap.parse_args()
    db_path = Path(a.db).resolve()
    if not db_path.exists():
        sys.exit(f"{db_path} not found")
    api = Api(a.url, a.user, Path(a.password_file).read_text(encoding="utf-8").strip())

    # read-only connection so analysis can never change curated data
    uri = f"sqlite:///file:{db_path.as_posix()}?mode=ro&uri=true"
    dbid = api.find("database", "database_name", DB_NAME)
    body = {"database_name": DB_NAME, "sqlalchemy_uri": uri, "expose_in_sqllab": True}
    if dbid:
        api.call("PUT", f"database/{dbid}", json=body)
    else:
        dbid = api.call("POST", "database/", json=body)["id"]

    dsid = api.find("dataset", "table_name", VIEW)
    if not dsid:
        dsid = api.call("POST", "dataset/", json={"database": dbid, "table_name": VIEW})["id"]
    api.call("PUT", f"dataset/{dsid}/refresh")

    ids = []
    for name, viz, params in CHARTS:
        p = {"viz_type": viz, "datasource": f"{dsid}__table", "adhoc_filters": [], "row_limit": 10000, **params}
        body = {"slice_name": name, "viz_type": viz, "datasource_id": dsid, "datasource_type": "table", "params": json.dumps(p)}
        cid = api.find("chart", "slice_name", name)
        if cid:
            api.call("PUT", f"chart/{cid}", json=body)
        else:
            cid = api.call("POST", "chart/", json=body)["id"]
        ids.append(cid)

    # dashboard layout: a 12-column grid, two charts per row (the big number on its own)
    pos = {"DASHBOARD_VERSION_KEY": "v2", "ROOT_ID": {"type": "ROOT", "id": "ROOT_ID", "children": ["GRID_ID"]},
           "GRID_ID": {"type": "GRID", "id": "GRID_ID", "children": [], "parents": ["ROOT_ID"]},
           "HEADER_ID": {"type": "HEADER", "id": "HEADER_ID", "meta": {"text": DASH}}}
    rows = [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]]
    for ri, row in enumerate(rows):
        rid = f"ROW-{ri}"
        pos["GRID_ID"]["children"].append(rid)
        pos[rid] = {"type": "ROW", "id": rid, "children": [], "parents": ["ROOT_ID", "GRID_ID"], "meta": {"background": "BACKGROUND_TRANSPARENT"}}
        for k in row:
            chid = f"CHART-{k}"
            pos[rid]["children"].append(chid)
            pos[chid] = {"type": "CHART", "id": chid, "children": [], "parents": ["ROOT_ID", "GRID_ID", rid],
                         "meta": {"chartId": ids[k], "width": 6, "height": 50, "sliceName": CHARTS[k][0]}}
    body = {"dashboard_title": DASH, "published": True, "position_json": json.dumps(pos)}
    did = api.find("dashboard", "dashboard_title", DASH)
    if did:
        api.call("PUT", f"dashboard/{did}", json=body)
    else:
        did = api.call("POST", "dashboard/", json=body)["id"]
    for cid in ids:  # attach charts to the dashboard
        api.call("PUT", f"chart/{cid}", json={"dashboards": [did]})
    print(f"database {dbid}, dataset {dsid}, charts {ids}, dashboard {did}")
    print(f"open: {a.url}/superset/dashboard/{did}/")


if __name__ == "__main__":
    main()
