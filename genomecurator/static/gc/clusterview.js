/* The shared view of the cluster: the same panel on the main server's page and on every sub-server's own page.
 *
 * A sub-server could show only the one job it was running, so from the machine actually doing the work there was no way to
 * see whether the analysis was moving or where the other chains had got to. Both pages now draw this one panel from what
 * `cluster.overview()` reports (the main server from its database, a sub-server through its own worker key).
 *
 *   renderCluster(element, overview)
 *
 * The styles are injected once here and written against BOTH pages' colour variables - the main pages inherit base.html's
 * (--line, --th, --stripe) and the sub-server page is a standalone document with its own (--bd, --bar, --card) - so neither
 * page has to learn the other's names. */
(function () {
  const S = `
.gcv{--gl:var(--bd,var(--line,#dde1e8));--gb:var(--bar,var(--th,#e6e9ef));--gc:var(--card,var(--bg2,#fff));--ga:var(--acc,#2b6cb0);--gok:var(--ok,#2f855a);--gbad:var(--bad,#c53030);--gwarn:var(--warn,#b7791f);--gmut:var(--mut,#6a7384)}
.gcv .gcv-cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(132px,1fr));gap:10px;margin:0 0 14px}
.gcv .gcv-card{border:1px solid var(--gl);border-radius:10px;padding:9px 12px;background:var(--gc)}
.gcv .gcv-card .k{font-size:12px;color:var(--gmut)}.gcv .gcv-card .v{font-size:19px;font-weight:700;line-height:1.25}
.gcv .gcv-card .s{font-size:11.5px;color:var(--gmut)}
.gcv .gcv-bar{height:7px;background:var(--gb);border-radius:4px;overflow:hidden;margin:5px 0 2px}
.gcv .gcv-bar i{display:block;height:100%;background:var(--ga);transition:width .4s}
.gcv .gcv-bar i.hi{background:var(--gwarn)}.gcv .gcv-bar i.crit{background:var(--gbad)}.gcv .gcv-bar i.job{background:var(--gok)}
.gcv .gcv-box{border:1px solid var(--gl);border-radius:10px;padding:12px 14px;margin:0 0 12px;background:var(--gc)}
.gcv .gcv-pipe{display:flex;align-items:center;flex-wrap:wrap;gap:6px;margin:8px 0 6px}
.gcv .gcv-step{border:1px solid currentColor;border-radius:14px;padding:2px 12px;font-size:12.5px;font-weight:600;color:var(--gmut)}
.gcv .gcv-step.done{color:var(--gok)}.gcv .gcv-step.now{color:var(--ga);box-shadow:0 0 0 3px color-mix(in srgb,var(--ga) 16%,transparent)}
.gcv .gcv-arrow{color:var(--gmut);font-size:12px}
.gcv .gcv-mach{display:grid;grid-template-columns:repeat(auto-fit,minmax(214px,1fr));gap:10px}
.gcv .gcv-m{border:1px solid var(--gl);border-radius:10px;padding:10px 12px;background:var(--gc)}
.gcv .gcv-m.you{border-color:var(--ga);box-shadow:0 0 0 2px color-mix(in srgb,var(--ga) 14%,transparent)}
.gcv .gcv-dot{display:inline-block;width:8px;height:8px;border-radius:50%;background:var(--gok);margin-right:6px;vertical-align:middle}
.gcv .gcv-dot.off{background:var(--gmut)}.gcv .gcv-dot.bad{background:var(--gbad)}
.gcv .gcv-j{border-top:1px solid var(--gl);padding:8px 2px}
.gcv .gcv-j:first-child{border-top:0}
.gcv .gcv-j.mine{background:color-mix(in srgb,var(--ga) 7%,transparent);border-radius:7px;padding-left:8px;padding-right:8px}
.gcv .gcv-mut{color:var(--gmut);font-size:12.5px}
.gcv h4{margin:16px 0 8px;font-size:14px}
.gcv .gcv-you{border:1px solid var(--ga);color:var(--ga);border-radius:10px;padding:0 7px;font-size:11px;font-weight:700;margin-left:6px;white-space:nowrap}
`;
  if (!document.getElementById("gcv-style")) {
    const t = document.createElement("style");
    t.id = "gcv-style"; t.textContent = S; document.head.appendChild(t);
  }
  const esc = s => String(s == null ? "" : s).replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const gb = b => b == null || b === 0 ? "-" : (b >= 2 ** 40 ? (b / 2 ** 40).toFixed(1) + " TB" : Math.round(b / 2 ** 30) + " GB");
  const dur = s => { s = Math.max(0, Math.round(s || 0)); const h = Math.floor(s / 3600), m = Math.floor(s % 3600 / 60), x = s % 60; return h ? h + "시간 " + m + "분" : m ? m + "분 " + x + "초" : x + "초"; };
  const ago = t => { if (!t) return "-"; const s = Math.round((Date.now() - new Date(t)) / 1000); return s < 90 ? s + "초 전" : s < 5400 ? Math.round(s / 60) + "분 전" : Math.round(s / 3600) + "시간 전"; };
  const bar = (pct, cls) => '<div class="gcv-bar"><i class="' + (cls || (pct >= 90 ? "crit" : pct >= 75 ? "hi" : "")) + '" style="width:' + Math.max(0, Math.min(100, pct || 0)) + '%"></i></div>';
  const STAGE = { bets: "시계 신호", compare: "모델 비교", dating: "본 분석", done: "완료", stopped: "중단" };

  function cards(T) {
    const c = [];
    const add = (k, v, s) => c.push('<div class="gcv-card"><div class="k">' + k + '</div><div class="v">' + v + '</div><div class="s">' + (s || "&nbsp;") + '</div></div>');
    add("서버", (T.nodes || 0) + "대", T.offline ? T.offline + "대 끊김" : "모두 연결됨");
    add("실행 중인 작업", (T.jobs || 0) + "개", T.threads_asked ? "요청 스레드 " + T.threads_asked : "");
    add("코어 사용", (T.cpus_used != null ? T.cpus_used : "-") + " / " + (T.cpus || "-"), T.cpu_pct != null ? T.cpu_pct + "%" : "");
    if (T.mem_total) add("메모리", gb(T.mem_used) + " / " + gb(T.mem_total), (T.mem_pct != null ? T.mem_pct + "% · " : "") + gb(T.mem_total - T.mem_used) + " 여유");
    add("디스크 여유", gb(T.disk_free), "");
    let warn = "";
    if (T.no_memory && T.no_memory.length) warn += '<div class="gcv-mut" style="color:var(--gwarn)">메모리를 보고하지 않는 서버: ' + esc(T.no_memory.join(", ")) + '</div>';
    if (T.other_code && T.other_code.length) warn += '<div class="gcv-mut" style="color:var(--gbad)">다른 코드를 쓰는 서버: ' + esc(T.other_code.join(", ")) + '</div>';
    if (T.no_build && T.no_build.length) warn += '<div class="gcv-mut" style="color:var(--gwarn)">코드를 알 수 없는 서버: ' + esc(T.no_build.join(", ")) + '</div>';
    return '<div class="gcv-cards">' + c.join("") + '</div>' + warn;
  }

  function planBox(p) {
    const g = p.progress, stages = (p.stages && p.stages.length ? p.stages : ["bets", "compare", "dating"]);
    const at = stages.indexOf(p.stage);
    const pipe = stages.map((st, i) => {
      const done = p.stage === "done" || (at >= 0 && i < at), now = p.stage === st;
      return '<span class="gcv-step ' + (done ? "done" : now ? "now" : "") + '">' + (STAGE[st] || st) + '</span>'
        + (i < stages.length - 1 ? '<span class="gcv-arrow">→</span>' : "");
    }).join("");
    let prog = "";
    if (g) {
      const pct = g.pct != null ? g.pct : (g.expected ? 100 * g.done / g.expected : 0);
      prog = bar(pct, "job") + '<div class="gcv-mut">' + (STAGE[p.stage] || p.stage) + ' <b>' + pct.toFixed(0) + '%</b> · 작업 <b>' + g.done + '/' + g.expected + '</b> 완료'
        + (g.running ? ' · 실행 중 ' + g.running : "") + (g.queued ? ' · 대기 ' + g.queued : "")
        + (g.failed ? ' · <span style="color:var(--gbad)">실패 ' + g.failed + '</span>' : "") + '</div>';
    }
    const ch = p.chosen ? '<div class="gcv-mut" style="margin-top:4px">선택된 모형 <b>' + esc(p.chosen.clock) + ' · ' + esc(p.chosen.tree_prior) + '</b>'
      + (p.chosen.runner_up ? ' (2위 ' + esc(p.chosen.runner_up) + ', 로그 BF ' + (p.chosen.log_bf_over_runner_up || 0).toFixed(1) + ')' : "") + '</div>' : "";
    return '<div class="gcv-box"><b>분석 계획 ' + esc(p.name) + '</b> <span class="gcv-mut">유전체 ' + (p.n_tips || "?") + '개 · ' + esc(p.model || "?")
      + ' · 격자 ' + (p.grid || 0) + '조합 · 사슬 ' + (p.chains || "?") + '개</span><div class="gcv-pipe">' + pipe + '</div>' + prog + ch
      + (p.log || []).map(l => '<div class="gcv-mut" style="font-size:12px">' + esc(l) + '</div>').join("") + '</div>';
  }

  function machine(n) {
    const m = n.metrics || {};
    const cpu = m.cpu_pct != null ? m.cpu_pct : (m.load1 != null && m.cpus ? 100 * m.load1 / m.cpus : null);
    const mem = m.mem_total ? 100 * m.mem_used / m.mem_total : null;
    return '<div class="gcv-m' + (n.you ? " you" : "") + '"><div><span class="gcv-dot ' + (n.online ? (n.paused ? "off" : "") : "bad") + '"></span><b>' + esc(n.name) + '</b>'
      + (n.you ? '<span class="gcv-you">이 서버</span>' : "") + '</div>'
      + '<div class="gcv-mut">' + (n.online ? (n.paused ? "정지됨" : "연결됨") : "끊김 · " + ago(n.last_seen)) + ' · 작업 ' + (n.jobs || 0) + '개'
      + (n.cpus ? ' · ' + n.cpus + '코어' + (n.cpus_free != null ? ' (여유 ' + n.cpus_free + ')' : "") : "") + '</div>'
      + (cpu != null ? bar(cpu) + '<div class="gcv-mut">' + (m.cpu_pct != null ? "CPU " + Math.round(cpu) + "%" : "부하 " + m.load1 + " / " + m.cpus + "코어") + '</div>' : "")
      + (mem != null ? bar(mem) + '<div class="gcv-mut">메모리 ' + gb(m.mem_used) + ' / ' + gb(m.mem_total) + '</div>' : '<div class="gcv-mut">메모리 보고 없음</div>')
      + '</div>';
  }

  function job(j) {
    const el = j.claimed ? (Date.now() - new Date(j.claimed)) / 1000 : 0;
    const eta = j.state === "running" && j.pct > 0.5 && el ? el * (100 - j.pct) / j.pct : null;
    return '<div class="gcv-j' + (j.mine ? " mine" : "") + '"><div style="display:flex;gap:8px;justify-content:space-between">'
      + '<span>' + esc(j.label) + (j.mine ? '<span class="gcv-you">이 서버</span>' : "") + '</span>'
      + '<b style="white-space:nowrap">' + (j.state === "running" ? (j.pct != null ? (+j.pct).toFixed(1) + "%" : "시작 중") : "대기") + '</b></div>'
      + (j.state === "running" ? bar(j.pct, "job") : "")
      + '<div class="gcv-mut">' + esc(j.node || "-") + (j.state === "running" ? ' · 수행 ' + dur(el) + ' · 남은 시간 ' + (eta ? "약 " + dur(eta)
        + " (완료 예상 " + new Date(Date.now() + eta * 1000).toLocaleString("ko-KR", { month: "numeric", day: "numeric", hour: "2-digit", minute: "2-digit" }) + ")" : "계산 중") : "")
      + (j.stale ? ' · <span style="color:var(--gbad)">응답 없음</span>' : "") + '</div></div>';
  }

  window.renderCluster = function (el, ov) {
    if (!el) return;
    if (!ov || ov.ok === false) {
      el.innerHTML = '<div class="gcv"><div class="gcv-box gcv-mut">' + esc((ov && ov.error) || "클러스터 현황을 읽지 못했습니다.") + '</div></div>';
      return;
    }
    const plans = ov.plans || [], jobs = ov.jobs || [];
    el.innerHTML = '<div class="gcv">' + cards(ov.totals || {})
      + (plans.length ? '<h4>분석 계획</h4>' + plans.map(planBox).join("") : "")
      + '<h4>서버</h4><div class="gcv-mach">' + (ov.nodes || []).map(machine).join("") + '</div>'
      + '<h4>실행 중인 작업 ' + jobs.filter(j => j.state === "running").length + '개'
      + (jobs.some(j => j.state === "queued") ? ' <span class="gcv-mut">· 대기 ' + jobs.filter(j => j.state === "queued").length + '개</span>' : "") + '</h4>'
      + '<div class="gcv-box" style="padding:4px 12px">' + (jobs.length ? jobs.map(job).join("") : '<div class="gcv-j gcv-mut">작업이 없습니다.</div>') + '</div>'
      + '</div>';
  };
})();
