// Venn diagram of the sources (NCBI, BV-BRC, Pathoplexus): how many genomes are in exactly which sources. Click a region to list its accessions.
(function () {
  const NAME = { ncbi: "NCBI", bvbrc: "BV-BRC", pathoplexus: "Pathoplexus" };
  const COL = { ncbi: "#3b82f6", bvbrc: "#f59e0b", pathoplexus: "#10b981" };
  const esc = s => String(s == null ? "" : s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const fmt = n => (n || 0).toLocaleString("ko-KR");
  // circle centres and the label position of every region, for 2 and 3 sources (classic, not area-proportional)
  const LAYOUT = {
    2: { c: [[120, 100], [200, 100]], r: 72, lab: { "0": [86, 100], "1": [234, 100], "0+1": [160, 100] }, name: [[80, 18], [240, 18]] },
    3: { c: [[125, 95], [205, 95], [165, 160]], r: 70, name: [[70, 22], [260, 22], [165, 245]],
         lab: { "0": [100, 80], "1": [230, 80], "2": [165, 195], "0+1": [165, 72], "0+2": [128, 140], "1+2": [202, 140], "0+1+2": [165, 115] } }
  };
  window.drawVenn = function (el, o, opt) {
    opt = opt || {};
    if (!o || !o.sources || !o.sources.length) { el.innerHTML = '<span class="mut">소스 기록이 없습니다.</span>'; return; }
    const S = o.sources, small = !!opt.small;
    if (S.length === 1) { el.innerHTML = `<div>${NAME[S[0]] || S[0]} <b>${fmt(o.sizes[S[0]])}</b>개 (다른 소스 없음)</div>`; return; }
    const L = LAYOUT[S.length], idx = s => S.indexOf(s);
    let g = "";
    S.forEach((s, i) => { const [x, y] = L.c[i]; g += `<circle cx="${x}" cy="${y}" r="${L.r}" fill="${COL[s]}" fill-opacity=".18" stroke="${COL[s]}" stroke-width="1.5"/>`; });
    S.forEach((s, i) => { const [x, y] = L.name[i]; g += `<text x="${x}" y="${y}" text-anchor="middle" font-size="12.5" font-weight="700" fill="${COL[s]}">${NAME[s] || s}</text><text x="${x}" y="${y + 14}" text-anchor="middle" font-size="11" fill="currentColor" fill-opacity=".75">${fmt(o.sizes[s])}개</text>`; });
    for (const r of o.regions) {
      const key = r.sources.map(idx).sort().join("+"), p = L.lab[key];
      if (!p) continue;
      const zero = !r.n, all = r.sources.length === S.length;
      g += `<g class="vr" data-key="${esc(r.key)}" style="cursor:${zero || !opt.onRegion ? "default" : "pointer"}"><title>${r.sources.map(s => NAME[s] || s).join(" ∩ ")}${all ? "" : " 에만"}: ${fmt(r.n)}개</title>` +
           `<circle cx="${p[0]}" cy="${p[1]}" r="17" fill="transparent"/><text x="${p[0]}" y="${p[1] + 4}" text-anchor="middle" font-size="${all ? 14 : 12.5}" font-weight="${all ? 700 : 600}" fill="currentColor" fill-opacity="${zero ? .35 : 1}">${fmt(r.n)}</text></g>`;
    }
    const H = S.length === 3 ? 262 : 182;
    const only = S.map(s => `${NAME[s] || s}에만 <b>${fmt(o.only[s])}</b>`).join(" · ");
    el.innerHTML = `<svg viewBox="0 0 330 ${H}" style="width:100%;max-width:${small ? 250 : 360}px;height:auto;display:block">${g}</svg>` +
      (small ? "" : `<div style="font-size:12.5px;margin-top:2px">합계 <b>${fmt(o.union)}</b>개 · 모든 소스에 <b>${fmt(o.in_all)}</b> · ${only}</div>`);
    if (opt.onRegion) el.querySelectorAll(".vr").forEach(n => n.addEventListener("click", () => {
      const r = o.regions.find(x => x.key === n.dataset.key); if (r && r.n) opt.onRegion(r);
    }));
  };
  window.vennRegionName = (r, nAll) => r.sources.map(s => NAME[s] || s).join(" ∩ ") + (r.sources.length < nAll ? "에만" : " 모두");
})();
