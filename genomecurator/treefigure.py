"""A publication-style figure of a BEAST summary tree, made by one function call so that it can be re-made from the same inputs.

Input: the annotated summary tree of a BEAST analysis (`beastprep.annotated_tree`: NEXUS with `posterior`, `height_95%_HPD`, and `country`, `lineage`, `date` on the
tips). Output: PNG, PDF and SVG plus a manifest (`repro.manifest`: input checksums, parameters, program / library / tool versions, checksums of the files written).

What the figure does, all of it set by the arguments and by nothing else:
- draws the tree on a calendar axis, children ordered so that the clade of interest (`focus`, a set of tip names) is last (a ladder), or by clade size;
- collapses side clades of `collapse_min` or more tips into triangles that list their composition by country; the clade of interest and the path to it are never collapsed;
- names every other tip with a label (`labels` maps a tip name to the text shown) and shows its country (colour), lineage (chip) and collection month;
- marks nodes with posterior support >= `support`, shades the clade of interest and draws the credible intervals given in `intervals` on it (MRCA) or on its stem.

Private genomes: a tip whose name starts with `private:`, or that is listed in `private_ids`, must have an alias in `labels`; otherwise the figure is refused
(the identifier of a private genome must not appear in a public output).

matplotlib is needed only here (`pip install genomecurator[figures]`).
"""
from __future__ import annotations

import collections
import datetime
import re
from pathlib import Path

from . import repro

LINEAGE_COLORS = ["#3b78b0", "#d1603d", "#5b9e55", "#8a5fa8", "#c9a227", "#2f9aa0"]
COUNTRY_COLORS = ["#8d8d8d", "#4daf4a", "#377eb8", "#984ea3", "#e41a1c", "#ff9f1a", "#a65628", "#17becf", "#bcbd22", "#f781bf", "#1b9e8f", "#aec7e8", "#7f7fff", "#ff7f7f"]
MIXED = "#8a8a8a"


# ---------------------------------------------------------------------------------------------------------------- the tree
def parse_nexus_tree(text: str) -> dict:
    """The first tree of a NEXUS file as nested dicts: {name, len, meta, children}. Handles a Translate table and [&key=value,...] annotations."""
    m = re.search(r"Translate\s*(.*?);\s*tree\s", text, re.S | re.I)
    tr = {}
    if m:
        for line in m.group(1).split(","):
            p = line.split(None, 1)
            if len(p) == 2:
                tr[p[0].strip()] = p[1].strip().strip("'\"")
    start = m.end() - 5 if m else 0
    s = re.search(r"tree\s+\S+\s*=\s*(?:\[&R\]\s*)?(.*);", text[start:], re.S).group(1)
    pos = [0]

    def meta():
        d = {}
        if pos[0] < len(s) and s[pos[0]] == "[":
            j = s.index("]", pos[0])
            body = s[pos[0] + 2:j]
            pos[0] = j + 1
            for key, val in re.findall(r'([A-Za-z0-9_.%]+)=(\{[^}]*\}|"[^"]*"|[^,]+)', body):
                d[key] = val.strip('"')
        return d

    def node():
        n = {"children": [], "meta": {}, "name": None, "len": 0.0}
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                n["children"].append(node())
                if s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                break
            pos[0] += 1
        else:
            j = pos[0]
            while s[j] not in "[:,)":
                j += 1
            n["name"] = tr.get(s[pos[0]:j], s[pos[0]:j])
            pos[0] = j
        n["meta"] = meta()
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            j = pos[0]
            while j < len(s) and s[j] not in ",)[":
                j += 1
            n["len"] = float(s[pos[0]:j])
            pos[0] = j
        n["meta"].update(meta())
        return n

    return node()


def tips(n: dict) -> list[dict]:
    return [n] if not n["children"] else [t for c in n["children"] for t in tips(c)]


def set_times(root: dict) -> None:
    """Calendar time (decimal year) of every node: the date of one dated tip minus its distance to the root, then branch lengths added down the tree."""
    def path(n, tgt, acc=0.0):
        if n is tgt:
            return acc
        for c in n["children"]:
            r = path(c, tgt, acc + c["len"])
            if r is not None:
                return r
    ref = next((t for t in tips(root) if t["meta"].get("date")), None)
    if ref is None:
        raise ValueError("The tree has no tip with a date; the summary tree must come from a dated analysis.")
    t0 = float(ref["meta"]["date"]) - path(root, ref)

    def down(n, t):
        n["t"] = t
        for c in n["children"]:
            down(c, t + c["len"])
    down(root, t0)


def mrca(n: dict, names: set) -> dict | None:
    if not names <= {t["name"] for t in tips(n)}:
        return None
    for c in n["children"]:
        r = mrca(c, names)
        if r:
            return r
    return n


def display_tree(root: dict, focus: set, collapse_min: int, order: str = "ladder") -> dict:
    """The tree as drawn: side clades of `collapse_min` or more tips become one 'tri' leaf; everything inside the clade of interest and on the path to it stays."""
    top = mrca(root, focus) if focus else None
    onpath = set()

    def mark(n):
        if top is not None and focus <= {t["name"] for t in tips(n)}:
            onpath.add(id(n))
            for c in n["children"]:
                mark(c)
    mark(root)
    inside = {id(x) for x in _all(top)} if top is not None else set()

    def build(n):
        ts = tips(n)
        if not n["children"]:
            return {"kind": "tip", "n": n, "t": n["t"]}
        if n is not root and id(n) not in onpath and id(n) not in inside and len(ts) >= collapse_min:
            return {"kind": "tri", "n": n, "t": n["t"], "tmax": max(x["t"] for x in ts), "tips": ts}
        kids = [build(c) for c in n["children"]]
        kids.sort(key=(lambda k: (id(k["n"]) in onpath, -len(tips(k["n"])))) if order == "ladder" else (lambda k: -len(tips(k["n"]))))
        return {"kind": "int", "n": n, "t": n["t"], "kids": kids}
    return build(root)


def _all(n):
    yield n
    for c in n["children"]:
        yield from _all(c)


def _leaves(d, out):
    if d["kind"] == "int":
        for k in d["kids"]:
            _leaves(k, out)
    else:
        out.append(d)
    return out


def ym(x: float) -> str:
    y = int(x)
    d = datetime.date(y, 1, 1) + datetime.timedelta(days=round((x - y) * 365.25))
    return f"{d.year}-{d.month:02d}"


def private_ids_of(run: dict) -> set:
    """The sequence names of the private genomes among the members of an analysis run (`private:<id>` in the member list, `<id>` in the FASTA and so in the tree)."""
    from . import datasets, nextstrain
    try:
        members = nextstrain.current_members(datasets.get(run["dataset_id"]))
    except Exception:
        members = nextstrain.run_members(run)
    return {m.split(":", 1)[1] for m in members if m.startswith("private:")}


def check_labels(tip_names, labels: dict, private_ids=()) -> None:
    """Refuse the figure when a private genome has no alias."""
    priv = set(private_ids)
    bad = [n for n in tip_names if (n.startswith("private:") or n in priv) and n not in labels]
    if bad:
        raise ValueError(f"{len(bad)} private genome(s) have no alias in `labels`; a private identifier must not appear in a figure (give each one an alias).")


# ---------------------------------------------------------------------------------------------------------------- the figure
def render(tree_file: str | Path, out_prefix: str | Path, *, focus: set | None = None, labels: dict | None = None, private_ids=(),
           intervals: list[dict] | None = None, title: str = "", subtitle: str = "", footnote: str = "", collapse_min: int = 3, order: str = "ladder",
           support: float = 0.95, formats=("png", "pdf", "svg"), width_in: float = 7.2, row_in: float = 0.147, dpi: int = 300, write_manifest: bool = True,
           extra_inputs=(), run: dict | None = None) -> dict:
    """Draw the figure. `run`: the analysis run the tree belongs to; its private genomes are then known and must have an alias in `labels`. `intervals`: [{"label", "lo", "hi", "color", "on": "clade" | "stem"}] (decimal years). -> {rows, collapsed, files, manifest}."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon

    focus, labels, intervals = set(focus or ()), dict(labels or {}), list(intervals or [])
    private_ids = set(private_ids) | (private_ids_of(run) if run else set())
    plt.rcParams.update({"font.family": ["Nimbus Sans", "Liberation Sans", "DejaVu Sans"], "font.size": 6.8, "axes.linewidth": 0.6, "svg.hashsalt": "genomecurator",
                         "svg.fonttype": "none", "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
    root = parse_nexus_tree(Path(tree_file).read_text(encoding="utf-8"))
    set_times(root)
    check_labels([t["name"] for t in tips(root)], labels, private_ids)
    D = display_tree(root, focus, collapse_min, order)
    rows = _leaves(D, [])
    for i, r in enumerate(rows):
        r["y"] = i

    def ypos(d):
        if d["kind"] != "int":
            return d["y"]
        d["y"] = sum(ypos(k) for k in d["kids"]) / len(d["kids"])
        return d["y"]
    ypos(D)
    lineages = sorted({t["meta"].get("lineage", "") for t in tips(root) if t["meta"].get("lineage")}, key=lambda x: (len(x), x))
    lc = {l: LINEAGE_COLORS[i % len(LINEAGE_COLORS)] for i, l in enumerate(lineages)}
    countries = sorted({t["meta"].get("country", "?") for t in tips(root)})
    cc = {c: COUNTRY_COLORS[i % len(COUNTRY_COLORS)] for i, c in enumerate(countries)}

    def branch_color(n):
        ls = {t["meta"].get("lineage") for t in tips(n)}
        return lc.get(next(iter(ls)), MIXED) if len(ls) == 1 else MIXED

    nr = len(rows)
    hh = row_in * nr + 3.85
    fig = plt.figure(figsize=(width_in, hh))
    ax = fig.add_subplot(fig.add_gridspec(1, 1, left=0.04, right=0.985, top=1 - 1.2 / hh, bottom=0.075)[0])
    t_all = [t["t"] for t in tips(root)]
    x0 = root["t"] - 0.15
    xl = max(t_all) + 0.45                                           # label column starts here
    ax.set_xlim(x0, xl + (max(t_all) - x0) * 0.95)
    ax.set_ylim(nr - 0.2, -2.2)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)

    def draw(d):
        for k in d["kids"]:
            col = branch_color(k["n"])
            ax.plot([d["t"]] * 2, [d["y"], k["y"]], color=col, lw=1.0, solid_capstyle="butt", zorder=2)
            ax.plot([d["t"], k["t"]], [k["y"]] * 2, color=col, lw=1.0, solid_capstyle="butt", zorder=2)
            if k["kind"] == "int":
                draw(k)
    draw(D)

    def post(d):
        if d["kind"] == "int":
            try:
                if float(d["n"]["meta"].get("posterior", 0)) >= support:
                    ax.scatter([d["t"]], [d["y"]], s=5, color="k", zorder=5, lw=0)
            except ValueError:
                pass
            for k in d["kids"]:
                post(k)
    post(D)
    focus_rows = [r for r in rows if r["kind"] == "tip" and r["n"]["name"] in focus]
    if focus_rows:
        ax.axhspan(min(r["y"] for r in focus_rows) - 0.5, max(r["y"] for r in focus_rows) + 0.5, color="#fbe3dc", zorder=0, lw=0)
    for r in rows:
        n = r["n"]
        if r["kind"] == "tip":
            m = n["meta"]
            hl = n["name"] in focus
            ax.scatter([r["t"]], [r["y"]], s=17 if hl else 11, color=cc.get(m.get("country", "?"), "#444"), edgecolor="k", linewidth=0.35, zorder=6)
            when = ym(float(m["date"])) if m.get("date") else ""
            txt = f"{labels.get(n['name'], n['name'])}   {m.get('country', '?')}   {when}"
            lin = m.get("lineage", "")
            ax.text(xl + 0.25, r["y"], txt, va="center", ha="left", fontsize=6.0, fontweight="bold" if hl else "normal", color="#b01010" if hl else "#222")
            ax.plot([r["t"] + 0.01, xl - 0.01], [r["y"]] * 2, color="#cfcfcf", lw=0.35, ls=(0, (1, 2)), zorder=1)
        else:
            ts = r["tips"]
            comp = collections.Counter(x["meta"].get("country", "?") for x in ts).most_common()
            lin = collections.Counter(x["meta"].get("lineage", "") for x in ts).most_common(1)[0][0]
            col = branch_color(n)
            ax.add_patch(Polygon([[r["t"], r["y"]], [r["tmax"], r["y"] - 0.42], [r["tmax"], r["y"] + 0.42]], closed=True, fc=col, ec="none", alpha=0.35, zorder=3))
            ax.plot([r["t"], r["tmax"]], [r["y"]] * 2, color=col, lw=0.6, zorder=4)
            ax.plot([r["tmax"] + 0.01, xl - 0.01], [r["y"]] * 2, color="#cfcfcf", lw=0.35, ls=(0, (1, 2)), zorder=1)
            ax.text(xl + 0.25, r["y"], f"{len(ts)} genomes: " + ", ".join(f"{k} {v}" for k, v in comp[:3]) + (", other" if len(comp) > 3 else ""),
                    va="center", ha="left", fontsize=5.8, color="#444", style="italic")
        ax.add_patch(plt.Rectangle((xl, r["y"] - 0.36), 0.2, 0.72, fc=lc.get(lin, MIXED), ec="none", zorder=4))
        ax.text(xl + 0.1, r["y"], lin, va="center", ha="center", fontsize=5.6, color="white", fontweight="bold", zorder=5)
    ax.text(xl + 0.1, -1.0, "lineage", fontsize=6, ha="center", color="#555")
    ax.text(xl + 0.25, -1.0, "genome   country   collection month", fontsize=6, ha="left", color="#555")

    if focus_rows:
        top = mrca(root, focus)

        def find(d, tgt, parent=None):
            if d["n"] is tgt:
                return d, parent
            if d["kind"] == "int":
                for k in d["kids"]:
                    r = find(k, tgt, d)
                    if r:
                        return r
        hit = find(D, top)
        if hit:
            node, parent = hit
            for iv in intervals:
                if iv.get("on", "clade") == "stem" and parent is not None:
                    ax.plot([iv["lo"], iv["hi"]], [parent["y"] + 0.62] * 2, color=iv.get("color", "#1f4e9c"), lw=2.4, solid_capstyle="butt", zorder=7)
                elif iv.get("on", "clade") == "clade":
                    ax.plot([iv["lo"], iv["hi"]], [node["y"] - 0.62] * 2, color=iv.get("color", "#c0182a"), lw=2.4, solid_capstyle="butt", zorder=7)

    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=cc[c], markeredgecolor="k", markeredgewidth=0.35, markersize=5, label=c) for c in countries]
    handles += [Line2D([0], [0], color=lc[l], lw=2.2, label=f"{l} branches") for l in lineages] + [Line2D([0], [0], color=MIXED, lw=2.2, label="mixed")]
    handles += [Line2D([0], [0], marker="<", color="w", markerfacecolor="#e8b6a6", markersize=7, label="collapsed clade")]
    if focus_rows:
        handles += [plt.Rectangle((0, 0), 1, 1, fc="#fbe3dc", ec="none", label="clade of interest")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=7, frameon=False, fontsize=6.0, columnspacing=0.9, handletextpad=0.5, handlelength=1.4)
    ticks = [yr + (mo - 1) / 12 for yr in range(int(x0), int(max(t_all)) + 2) for mo in (1, 4, 7, 10) if x0 <= yr + (mo - 1) / 12 <= max(t_all) + 0.1]
    ax.set_xticks(ticks)
    ax.set_xticklabels([ym(t + 0.001) for t in ticks], fontsize=6.4)
    if ticks:
        ax.spines["bottom"].set_bounds(ticks[0], ticks[-1])
    if title:
        fig.text(0.04, 1 - 0.13 / hh, title, fontsize=8.6, fontweight="bold", va="top")
    if subtitle:
        fig.text(0.04, 1 - 0.34 / hh, subtitle, fontsize=6.0, color="#333", va="top", linespacing=1.4)
    if footnote:
        fig.text(0.04, 0.014, footnote, fontsize=5.4, color="#555", va="bottom", linespacing=1.4)

    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    files = []
    for fmt in formats:
        f = out_prefix.with_suffix("." + fmt)
        meta = {"png": {"Software": None}, "pdf": {"CreationDate": None, "Creator": "genomecurator"}, "svg": {"Date": None}}.get(fmt, {})
        fig.savefig(f, dpi=dpi, metadata=meta)
        files.append(str(f))
    plt.close(fig)
    params = {"focus": sorted(focus), "collapse_min": collapse_min, "order": order, "support": support, "width_in": width_in, "row_in": row_in, "dpi": dpi,
              "intervals": intervals, "title": title, "subtitle": subtitle, "footnote": footnote, "labels_for": sorted(labels), "rows": nr}
    man = repro.manifest(files, [tree_file, *extra_inputs], params, write=write_manifest)
    return {"rows": nr, "collapsed": sum(1 for r in rows if r["kind"] == "tri"), "files": files, "manifest": man}
