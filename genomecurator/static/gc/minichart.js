// Small history charts for the Cluster pages (no library): CPU, memory, disk and network from a Netdata history, with job periods shaded.
(function () {
  const esc = s => String(s == null ? "" : s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const hm = t => { const d = new Date(t * 1000); return String(d.getHours()).padStart(2, "0") + ":" + String(d.getMinutes()).padStart(2, "0"); };
  function one(t, series, unit, fixedMax, jobs) {
    const W = 600, H = 120, L = 34, R = 6, T = 8, B = 18;
    const vals = series.flatMap(s => (s.v || []).filter(x => x != null));
    if (!t.length || !vals.length) return '<div class="mut" style="font-size:12px">자료 없음</div>';
    const max = fixedMax || Math.max(1e-9, ...vals) * 1.15, t0 = t[0], t1 = t[t.length - 1] || t0 + 1;
    const X = x => L + (W - L - R) * (x - t0) / Math.max(1, t1 - t0), Y = y => T + (H - T - B) * (1 - Math.min(y, max) / max);
    let g = "";
    for (const j of jobs || []) {
      const a = Math.max(t0, j.start), b = Math.min(t1, j.end || t1);
      if (b > a) g += `<rect x="${X(a)}" y="${T}" width="${Math.max(1, X(b) - X(a))}" height="${H - T - B}" fill="var(--jobband,rgba(72,187,120,.16))"><title>${esc(j.label)}</title></rect>`;
    }
    for (const f of [0, 0.5, 1]) g += `<line x1="${L}" x2="${W - R}" y1="${Y(max * f)}" y2="${Y(max * f)}" stroke="currentColor" stroke-opacity=".12"/><text x="${L - 4}" y="${Y(max * f) + 3}" text-anchor="end" font-size="9" fill="currentColor" fill-opacity=".6">${+(max * f).toFixed(max * f < 10 ? 1 : 0)}</text>`;
    for (const s of series) {
      const pts = t.map((x, i) => s.v && s.v[i] != null ? `${X(x).toFixed(1)},${Y(s.v[i]).toFixed(1)}` : null).filter(Boolean).join(" ");
      g += `<polyline points="${pts}" fill="none" stroke="${s.c}" stroke-width="1.6"/>`;
    }
    g += `<text x="${L}" y="${H - 4}" font-size="9" fill="currentColor" fill-opacity=".6">${hm(t0)}</text><text x="${W - R}" y="${H - 4}" text-anchor="end" font-size="9" fill="currentColor" fill-opacity=".6">${hm(t1)}</text>`;
    const last = series.map(s => { const v = (s.v || []).filter(x => x != null).pop(); return `<span style="color:${s.c}">■</span> ${esc(s.n)} ${v == null ? "-" : v}${unit}`; }).join(" &nbsp;");
    return `<svg viewBox="0 0 ${W} ${H}" style="width:100%;height:auto;display:block">${g}</svg><div style="font-size:12px;opacity:.8">${last}</div>`;
  }
  window.miniCharts = function (el, h) {
    if (!h || !h.ok) { el.innerHTML = '<div class="mut">' + esc(h && h.error ? "Netdata에서 기록을 읽지 못했습니다: " + h.error : "이 서버에 Netdata가 없어 기록 그래프가 없습니다(현재 값만 표시).") + "</div>"; return; }
    const jobs = h.jobs || [], box = (title, inner) => `<div><div style="font-size:12.5px;font-weight:600;margin-bottom:2px">${title}</div>${inner}</div>`;
    el.innerHTML = `<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px">` +
      box("CPU (%)", one(h.t, [{ n: "사용", v: h.cpu, c: "#3b82f6" }], "%", 100, jobs)) +
      box("메모리 (%)", one(h.t, [{ n: "사용", v: h.mem, c: "#a855f7" }], "%", 100, jobs)) +
      box("디스크 (MB/s)", one(h.t, [{ n: "읽기", v: h.disk_read, c: "#10b981" }, { n: "쓰기", v: h.disk_write, c: "#f59e0b" }], "", 0, jobs)) +
      box("네트워크 (Mbit/s)", one(h.t, [{ n: "받기", v: h.net_in, c: "#06b6d4" }, { n: "보내기", v: h.net_out, c: "#ef4444" }], "", 0, jobs)) +
      `</div><div class="mut" style="font-size:12px;margin-top:4px">초록 띠: 이 서버에서 클러스터 작업이 돈 구간 · 출처 ${esc(h.source || "")}</div>`;
  };
})();
