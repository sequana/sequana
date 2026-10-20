// Multi-resolution (min/max/mean) genome-track viewer: single page, one sequence.
// All data comes from the JSON in <script id="meta"> (see genome_track/builder.py).
(function () {
  "use strict";
  const META = JSON.parse(document.getElementById("meta").textContent);
  const $ = (id) => document.getElementById(id);
  const status = (msg, err) => { const s = $("status"); s.textContent = msg; s.className = err ? "err" : ""; };
  window.addEventListener("error", (e) => status("Error: " + e.message, true));
  window.addEventListener("unhandledrejection", (e) => status("Error: " + (e.reason && e.reason.message || e.reason), true));

  const POS0 = META.pos0, N = META.length, POS1 = POS0 + N - 1;
  const TH = META.thresholds; // optional confidence band: (threshold * sigma + mu) * baseline
  const LBL = META.labels;
  const COLORS = { low: "rgba(0,114,178,0.28)", high: "rgba(213,94,0,0.28)", lowEdge: "#0072b2", highEdge: "#d55e00" };
  const fmt = (v) => Math.round(v).toLocaleString("en-US");
  const pct = (u, v) => v.map((t) => Math.round(t * 100) + "%");
  const AXIS = 60; // same left/right axis width on every plot keeps their x axes aligned

  // ---- decoding -----------------------------------------------------------
  async function inflate(b64) {
    const bin = atob(b64);
    const bytes = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("deflate-raw"));
    return new Uint8Array(await new Response(stream).arrayBuffer());
  }
  function unshuffle(u8, n) { // 4 byte planes -> float32
    const out = new Uint8Array(n * 4);
    for (let b = 0; b < 4; b++) {
      const plane = u8.subarray(b * n, (b + 1) * n);
      for (let i = 0; i < n; i++) out[i * 4 + b] = plane[i];
    }
    return new Float32Array(out.buffer);
  }
  async function decodeLevels() {
    const levels = [];
    for (const L of META.levels) {
      const level = { bin: L.bin, n: L.n };
      for (const [name, b64] of Object.entries(L.arrays)) level[name] = unshuffle(await inflate(b64), L.n);
      if (L.bin === 1) { level.mean = level.min = level.max = level.cov; }
      levels.push(level);
    }
    GCL = levels.find((l) => l.gc) || null;
    return levels;
  }

  // ---- view helpers -------------------------------------------------------
  let LEVELS = [], GCL = null; // GCL: finest level that stores the GC content
  function chooseLevel(x0, x1, widthPx) {
    for (const L of LEVELS) if ((x1 - x0) / L.bin <= 3 * widthPx) return L;
    return LEVELS[LEVELS.length - 1];
  }
  const thrHigh = (rm) => (TH.high * TH.sigma + TH.mu) * rm;
  const thrLow = (rm) => (TH.low * TH.sigma + TH.mu) * rm;
  const ICOL_GC = TH ? 7 : 5; // column of the GC series in the main plot data

  function slice(L, x0, x1) {
    const i0 = Math.max(0, Math.floor((x0 - POS0) / L.bin) - 1);
    const i1 = Math.min(L.n, Math.ceil((x1 - POS0) / L.bin) + 1);
    const n = Math.max(0, i1 - i0);
    const x = new Float64Array(n), hi = TH ? new Float64Array(n) : null, lo = TH ? new Float64Array(n) : null;
    const rm = L.rm.subarray(i0, i1);
    let ymax = 0;
    for (let k = 0; k < n; k++) {
      x[k] = POS0 + (i0 + k) * L.bin + (L.bin - 1) / 2;
      if (TH) { hi[k] = thrHigh(rm[k]); lo[k] = thrLow(rm[k]); }
    }
    const mx = L.max.subarray(i0, i1), mn = L.min.subarray(i0, i1);
    let ymin = Infinity, thrMax = 0, thrMin = Infinity;
    for (let k = 0; k < n; k++) {
      if (mx[k] > ymax) ymax = mx[k];
      if (mn[k] < ymin) ymin = mn[k];
      if (TH) { if (hi[k] > thrMax) thrMax = hi[k]; if (lo[k] < thrMin) thrMin = lo[k]; }
    }
    let gc = null;
    if (GCL) {
      if (L.gc) gc = L.gc.subarray(i0, i1);
      else { // finer than the stored GC resolution: reuse the nearest stored level
        gc = new Float32Array(n);
        for (let k = 0; k < n; k++) gc[k] = GCL.gc[Math.min(GCL.n - 1, Math.floor(((i0 + k) * L.bin) / GCL.bin))];
      }
    }
    const data = [x, L.mean.subarray(i0, i1), mn, mx, rm];
    if (TH) data.push(hi, lo);
    if (gc) data.push(gc);
    return { data, ymax, ymin, thrMax, thrMin, level: L, gc };
  }

  // ---- ROIs ---------------------------------------------------------------
  // [start, end, size, mean_cov, mean_rm, log2_ratio, max_zscore, type]
  const ROIS = META.rois.map((r) => ({ start: r[0], end: r[1], size: r[2], mean_cov: r[3], mean_rm: r[4], log2: r[5], maxz: r[6], type: r[7], genes: r[8] || [] }));
  let selected = null;

  // ---- annotation (optional) ----------------------------------------------
  // [start, end, strand, lane, type, name, product, id]  (1-based, end inclusive)
  const ANN = META.ann ? META.ann.features.map((f) => ({ start: f[0], end: f[1], strand: f[2], lane: f[3], type: f[4], name: f[5], product: f[6], id: f[7] })) : null;
  const LANE_H = 16;
  let selectedFeature = null;
  const geneNames = (r) => r.genes.map((i) => ANN[i].name);
  const geneText = (r, max) => {
    const n = geneNames(r);
    return n.length ? (n.length > max ? n.slice(0, max).join(", ") + ` +${n.length - max}` : n.join(", ")) : "";
  };
  const productText = (r) => [...new Set(r.genes.map((i) => ANN[i].product).filter(Boolean))].slice(0, 2).join("; ");
  if (ANN) for (const r of ROIS) { r.geneStr = geneText(r, 99); r.geneShort = geneText(r, 3); r.productStr = productText(r); }

  function drawRois(u, ids) {
    const { ctx } = u, { top, height, left, width } = u.bbox;
    const xmin = u.scales.x.min, xmax = u.scales.x.max;
    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, width, height); ctx.clip();
    const dpr = window.devicePixelRatio || 1;
    for (const r of ROIS) {
      if (r.end < xmin || r.start > xmax) continue;
      const a = u.valToPos(r.start, "x", true), b = u.valToPos(r.end, "x", true);
      if (b - a >= 3 * dpr) { // wide enough: shade the full height
        ctx.fillStyle = COLORS[r.type];
        ctx.fillRect(a, top, b - a, height);
      } else { // narrow: a small tick at the bottom keeps the coverage readable
        ctx.fillStyle = r.type === "low" ? COLORS.lowEdge : COLORS.highEdge;
        ctx.fillRect(a - dpr / 2, top + height - 9 * dpr, Math.max(b - a, 1.5 * dpr), 9 * dpr);
      }
      if (r === selected) {
        ctx.strokeStyle = "#111"; ctx.lineWidth = 2 * dpr;
        ctx.strokeRect(a, top, Math.max(b - a, 3 * dpr), height);
      }
    }
    ctx.restore();
  }

  // Share only the crosshair between the coverage and the GC plots: by default uPlot also
  // mirrors drag selections, which would zoom the other plot a second time.
  const onlyCrosshair = (type) => type === "mousemove" || type === "mouseleave";
  const CROSSHAIR_SYNC = { key: "cov", filters: { pub: onlyCrosshair, sub: onlyCrosshair } };

  function drawAnnotation(u) {
    const { ctx } = u, { left, top, width, height } = u.bbox, dpr = window.devicePixelRatio || 1;
    const xmin = u.scales.x.min, xmax = u.scales.x.max;
    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, width, height); ctx.clip();
    ctx.font = `${11 * dpr}px system-ui, sans-serif`; ctx.textBaseline = "middle";
    for (const f of ANN) {
      if (f.end < xmin || f.start > xmax) continue;
      const a = u.valToPos(f.start - 0.5, "x", true), b = u.valToPos(f.end + 0.5, "x", true);
      const w = Math.max(b - a, dpr), y = top + (f.lane * LANE_H + 2) * dpr, h = (LANE_H - 4) * dpr;
      ctx.fillStyle = f.type === "gene" ? (f.strand === "-" ? "#8a63a8" : "#3f6fa6") : "#7a7a7a";
      ctx.fillRect(a, y, w, h);
      if (w > 12 * dpr) { // strand arrow at the downstream end
        const dir = f.strand === "-" ? -1 : f.strand === "+" ? 1 : 0, tip = dir > 0 ? a + w : a;
        if (dir) {
          ctx.fillStyle = "rgba(255,255,255,0.85)"; ctx.beginPath();
          ctx.moveTo(tip - dir * 1 * dpr, y + h / 2); ctx.lineTo(tip - dir * 6 * dpr, y + 2 * dpr); ctx.lineTo(tip - dir * 6 * dpr, y + h - 2 * dpr); ctx.fill();
        }
      }
      if (w > 56 * dpr) { // label when there is room
        ctx.save(); ctx.beginPath(); ctx.rect(a, y, w - 8 * dpr, h); ctx.clip();
        ctx.fillStyle = "#fff"; ctx.fillText(f.name, a + 3 * dpr, y + h / 2 + dpr); ctx.restore();
      }
      if (f === selectedFeature) { ctx.strokeStyle = "#111"; ctx.lineWidth = 2 * dpr; ctx.strokeRect(a, y, w, h); }
    }
    // the selected region of interest, if any
    if (selected && selected.end >= xmin && selected.start <= xmax) {
      const a = u.valToPos(selected.start, "x", true), b = u.valToPos(selected.end, "x", true);
      ctx.fillStyle = selected.type === "low" ? COLORS.low : COLORS.high;
      ctx.fillRect(a, top, Math.max(b - a, 2 * dpr), height);
    }
    ctx.restore();
  }

  // ---- plots --------------------------------------------------------------
  let main, gcPlot = null, annPlot = null, busy = false, ymaxCap = Infinity, hashTimer = null;
  const widthOf = () => Math.max(300, $("plot").clientWidth);

  function clipOn() { return $("clip").checked; }
  function yTop(ymax) { return Math.max(5, clipOn() ? Math.min(ymax * 1.05, ymaxCap) : ymax * 1.05); }
  // Whole chromosome: 0 .. normal range. Zoomed: fit the visible coverage and thresholds.
  function yRange(s, min, max) {
    const zoomed = $("fity").checked && (max - min) < 0.8 * N;
    if (!zoomed) return [0, yTop(s.ymax)];
    let lo = Math.min(s.ymin, s.thrMin), hi = Math.max(s.ymax, s.thrMax);
    if (clipOn()) hi = Math.min(hi, Math.max(ymaxCap, s.ymax));
    if (!(hi > lo)) { lo = Math.max(0, lo - 1); hi = lo + 2; }
    const pad = 0.06 * (hi - lo);
    return [Math.max(0, lo - pad), hi + pad];
  }

  function refreshMain() {
    const { min, max } = main.scales.x;
    if (!(max > min)) return; // scales not committed yet (a refresh is scheduled by their setScale hook)
    drawn.min = min; drawn.max = max;
    const L = chooseLevel(min, max, main.bbox.width / (window.devicePixelRatio || 1));
    const s = slice(L, min, max);
    busy = true;
    main.setData(s.data, true); // true: re-applies the current x range, so the drawn index window follows the new data
    const [ya, yb] = manualY || yRange(s, min, max);
    main.setScale("y", { min: ya, max: yb });
    if (META.has_gc) main.setScale("gc", { min: 0, max: 1 }); // right axis: always 0 - 100 %
    if (gcPlot && $("gcpanel").checked) {
      gcPlot.setData([s.data[0], s.gc], true); // first the data (it re-applies the old x range) ...
      gcPlot.setScale("x", { min, max });       // ... then the new range
      gcPlot.setScale("y", { min: 0, max: 1 }); // always the full 0 - 100 % range
    }
    if (annPlot) annPlot.setScale("x", { min, max });
    busy = false;
    status(`view ${fmt(min)} - ${fmt(max)} bp (${fmt(max - min)} bp), resolution: ${L.bin === 1 ? "1 base" : L.bin + " bp bins (min/max/mean)"}, ${s.data[0].length.toLocaleString("en-US")} points drawn`);
    clearTimeout(hashTimer); // browsers throttle history updates: write the hash once the view settles
    hashTimer = setTimeout(() => history.replaceState(null, "", `#x=${Math.round(min)}-${Math.round(max)}`), 200);
  }
  let refreshQueued = false;
  const drawn = { min: NaN, max: NaN }; // x range of the data currently in the main plot
  function scheduleRefresh() {
    if (refreshQueued) return;
    refreshQueued = true;
    queueMicrotask(() => { refreshQueued = false; refreshMain(); });
  }
  // y range chosen by a rectangle selection; any other navigation goes back to the fitted y axis
  let manualY = null;
  function setRange(min, max, keepY) {
    if (!keepY) manualY = null;
    min = Math.max(POS0, min); max = Math.min(POS1, max);
    if (max - min < 20) { const c = (min + max) / 2; min = c - 10; max = c + 10; }
    main.setScale("x", { min, max });
  }

  // mouse wheel: zoom around the pointer (works on every plot sharing the x axis).
  // One notch (deltaY = 100) changes the width by about 30%. Events are accumulated
  // and applied once per animation frame.
  let wheelTarget = null, wheelFrame = 0;
  function attachWheel(u) {
    u.over.addEventListener("wheel", (e) => {
      e.preventDefault();
      const left = e.clientX - u.over.getBoundingClientRect().left;
      const at = u.posToVal(left, "x");
      const cur = wheelTarget || { min: main.scales.x.min, max: main.scales.x.max };
      const delta = e.deltaY * (e.deltaMode === 1 ? 33 : 1); // lines -> pixels (Firefox)
      const factor = Math.exp(delta * 0.0035); // < 1 zooms in, > 1 zooms out
      let a = at - (at - cur.min) * factor, b = at + (cur.max - at) * factor;
      if (b - a >= POS1 - POS0) { a = POS0; b = POS1; }
      else if (a < POS0) { b += POS0 - a; a = POS0; }
      else if (b > POS1) { a -= b - POS1; b = POS1; }
      wheelTarget = { min: a, max: b };
      if (!wheelFrame) wheelFrame = requestAnimationFrame(() => {
        const t = wheelTarget; wheelTarget = null; wheelFrame = 0;
        setRange(t.min, t.max);
      });
    }, { passive: false });
  }

  // Panning: drag with the right mouse button (or Shift + left button, e.g. on a trackpad).
  // A capture listener on the plot root keeps uPlot from starting its own rectangle selection.
  let panFrame = 0, panTarget = null;
  function attachPan(u) {
    const over = u.over;
    over.addEventListener("contextmenu", (e) => e.preventDefault()); // the right button pans
    u.root.addEventListener("mousedown", (e) => {
      const wanted = e.button === 2 || (e.button === 0 && e.shiftKey);
      if (!wanted || !over.contains(e.target)) return;
      e.preventDefault(); e.stopImmediatePropagation(); // no uPlot selection
      const rect = over.getBoundingClientRect();
      const start = { x: e.clientX, y: e.clientY, min: main.scales.x.min, max: main.scales.x.max, ymin: main.scales.y.min, ymax: main.scales.y.max, yFixed: !!manualY };
      over.classList.add("panning");
      const move = (ev) => {
        const span = start.max - start.min;
        let a = start.min - ((ev.clientX - start.x) / rect.width) * span, b = a + span;
        if (a < POS0) { a = POS0; b = a + span; } else if (b > POS1) { b = POS1; a = b - span; }
        let y = null;
        if (start.yFixed) { // a rectangle fixed the y range: pan vertically too
          const dy = ((ev.clientY - start.y) / rect.height) * (start.ymax - start.ymin);
          y = [start.ymin + dy, start.ymax + dy];
        }
        panTarget = { a, b, y };
        if (!panFrame) panFrame = requestAnimationFrame(() => {
          const t = panTarget; panFrame = 0;
          if (t.y) manualY = t.y;
          setRange(t.a, t.b, !!t.y);
        });
      };
      const up = () => { over.classList.remove("panning"); window.removeEventListener("mousemove", move); window.removeEventListener("mouseup", up); };
      window.addEventListener("mousemove", move);
      window.addEventListener("mouseup", up);
    }, true);
  }

  // tooltip on hover, click to zoom to a feature
  function featureAt(u, clientX, clientY) {
    const r = u.over.getBoundingClientRect();
    const bp = u.posToVal(clientX - r.left, "x"), tol = (3 * (u.scales.x.max - u.scales.x.min)) / r.width;
    const lane = Math.floor((clientY - r.top) / LANE_H);
    let best = null;
    for (const f of ANN) if (f.lane === lane && f.start - tol <= bp && f.end + tol >= bp) { if (!best || f.end - f.start < best.end - best.start) best = f; }
    return best;
  }
  function wireAnnotation(u) {
    const tip = $("tip"); let down = null;
    const esc = (t) => String(t).replace(/[&<>]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));
    u.over.addEventListener("mousemove", (e) => {
      const f = featureAt(u, e.clientX, e.clientY);
      u.over.style.cursor = f ? "pointer" : "";
      if (!f) { tip.style.display = "none"; return; }
      tip.innerHTML = `<b>${esc(f.name)}</b> <span style="opacity:.7">${esc(f.type)} (${f.strand})</span><br>${esc(f.product || "no product")}<br>` +
        `${fmt(f.start)} - ${fmt(f.end)} (${fmt(f.end - f.start + 1)} bp)<br><span style="opacity:.7">${esc(f.id)}</span>`;
      tip.style.display = "block";
      tip.style.left = Math.min(e.pageX + 14, document.documentElement.clientWidth - 330) + "px"; tip.style.top = e.pageY + 14 + "px";
    });
    u.over.addEventListener("mouseleave", () => { tip.style.display = "none"; });
    u.over.addEventListener("mousedown", (e) => { down = [e.clientX, e.clientY]; });
    u.over.addEventListener("click", (e) => { // a click (not a drag) on a feature zooms to it
      if (!down || Math.abs(e.clientX - down[0]) > 3 || Math.abs(e.clientY - down[1]) > 3) return;
      const f = featureAt(u, e.clientX, e.clientY);
      if (!f) return;
      selectedFeature = f;
      const pad = Math.max(500, (f.end - f.start) * 0.5);
      setRange(f.start - pad, f.end + pad);
      u.redraw(false);
    });
  }

  function build() {
    const fullRange = [POS0, POS1];
    let x0 = fullRange[0], x1 = fullRange[1];
    const m = /#x=(\d+)-(\d+)/.exec(location.hash);
    if (m) { x0 = +m[1]; x1 = +m[2]; }

    // y cap used when "clip" is on: 2 x the highest threshold (1.5 x the baseline without thresholds), or 10 x the mean
    const L0 = LEVELS[LEVELS.length - 1];
    let maxHi = 0; for (let i = 0; i < L0.n; i++) maxHi = Math.max(maxHi, TH ? thrHigh(L0.rm[i]) : 1.5 * L0.rm[i]);
    ymaxCap = Math.min(2 * maxHi, 10 * META.signal_mean) || Infinity;

    const w = widthOf();
    const s = slice(chooseLevel(x0, x1, w), x0, x1);
    const mkSeries = (label, stroke, extra) => Object.assign({ label, stroke, width: 1, points: { show: false } }, extra || {});
    main = new uPlot({
      width: w, height: 420,
      // y is not auto-ranged: refreshMain() sets it explicitly (zoomed vs whole chromosome)
      scales: Object.assign(
        { x: { time: false, auto: false, range: (u, a, b) => [a, b] }, y: { auto: false, range: (u, a, b) => [a, b] } },
        META.has_gc ? { gc: { auto: false, range: (u, a, b) => [a, b] } } : {}
      ),
      axes: [{ label: "Position (bp)", values: (u, v) => v.map((t) => fmt(t)) }, { label: LBL.y_axis, size: AXIS }].concat(
        META.has_gc ? [{ scale: "gc", side: 1, label: LBL.gc, size: AXIS, grid: { show: false }, stroke: "#007a5a", values: pct }] : []
      ),
      series: [
        { label: "Position" },
        mkSeries(LBL.signal, "#111"),
        mkSeries("min", "rgba(0,0,0,0)", { width: 0 }),
        mkSeries("max", "rgba(0,0,0,0)", { width: 0 }),
        mkSeries(LBL.baseline, "#0099cc"),
      ].concat(TH ? [mkSeries("Thresholds (high)", "#d62728", { dash: [6, 4] }), mkSeries("Thresholds (low)", "#d62728", { dash: [6, 4] })] : [])
        .concat(META.has_gc ? [mkSeries(LBL.gc, "rgba(0,158,115,0.85)", { scale: "gc", width: 1.5 })] : []),
      bands: [{ series: [3, 2], fill: "rgba(60,60,60,0.22)" }],
      legend: { mount: (u, el) => $("legend").appendChild(el) }, // live values above the plot, genes right under it
      cursor: {
        drag: { x: true, y: true, setScale: false }, // rectangle selection, handled in the setSelect hook
        sync: META.has_gc ? CROSSHAIR_SYNC : undefined,
        bind: { dblclick: (u, targ, handler) => (e) => { setRange(POS0, POS1); return null; } },
      },
      hooks: {
        drawClear: [(u) => drawRois(u)],
        setSelect: [(u) => { // rectangle: zoom on x and y; a (nearly) flat drag zooms on x only
          const { left, top, width, height } = u.select;
          if (width < 3) return;
          const a = u.posToVal(left, "x"), b = u.posToVal(left + width, "x");
          const yhi = u.posToVal(top, "y"), ylo = u.posToVal(top + height, "y");
          u.setSelect({ left: 0, top: 0, width: 0, height: 0 }, false);
          if (height >= 8) {
            const pad = Math.max(1, (yhi - ylo) * 0.02);
            manualY = [Math.max(0, ylo - pad), yhi + pad];
            setRange(a, b, true);
            scheduleRefresh(); // after the x scale is committed; also applies y when x is unchanged
          } else setRange(a, b);
        }],
        // do not change scales from inside the hook (uPlot is still committing them): defer
        setScale: [(u, key) => {
          // skip when the x range is the one already drawn (setData() re-applies the scale: no loop)
          if (key === "x" && !busy && (u.scales.x.min !== drawn.min || u.scales.x.max !== drawn.max)) scheduleRefresh();
        }],
      },
    }, s.data, $("plot"));
    main.setScale("x", { min: x0, max: x1 }); // the setScale hook schedules the first refresh

    // GC track: own panel below the coverage, x axis kept in sync with the main plot
    if (META.has_gc) {
      gcPlot = new uPlot({
        width: w, height: 130,
        scales: { x: { time: false, auto: false, range: (u, a, b) => [a, b] }, y: { auto: false, range: (u, a, b) => [a, b] } },
        axes: [{ values: (u, v) => v.map((t) => fmt(t)) }, { label: LBL.gc, size: AXIS, values: pct }, { side: 1, size: AXIS, label: " ", values: () => [], grid: { show: false }, ticks: { show: false } }],
        cursor: { drag: { x: true, y: false, setScale: false }, sync: CROSSHAIR_SYNC },
        series: [{ label: "Position" }, mkSeries(LBL.gc, "#009e73", { width: 1.5 })],
        hooks: {
          setSelect: [(u) => {
            if (u.select.width < 3) return;
            const a = u.posToVal(u.select.left, "x"), b = u.posToVal(u.select.left + u.select.width, "x");
            u.setSelect({ left: 0, top: 0, width: 0, height: 0 }, false);
            setRange(a, b);
          }],
        },
      }, [s.data[0], s.gc], $("gcplot"));
      attachWheel(gcPlot);
    }
    if (ANN) {
      const lanes = META.ann.lanes;
      annPlot = new uPlot({
        width: w, height: lanes * LANE_H + 8,
        padding: [0, null, 0, null],
        legend: { show: false },
        scales: { x: { time: false, auto: false, range: (u, a, b) => [a, b] }, y: { auto: false, range: [0, 1] } },
        axes: [{ show: false }, { label: "Genes", size: AXIS, values: () => [], grid: { show: false }, ticks: { show: false } },
               { side: 1, size: AXIS, label: " ", values: () => [], grid: { show: false }, ticks: { show: false } }],
        series: [{}, { stroke: "rgba(0,0,0,0)", points: { show: false } }],
        cursor: { drag: { x: true, y: false, setScale: false }, sync: CROSSHAIR_SYNC, points: { show: false } },
        hooks: {
          draw: [drawAnnotation],
          setSelect: [(u) => {
            if (u.select.width < 3) return;
            const a = u.posToVal(u.select.left, "x"), b = u.posToVal(u.select.left + u.select.width, "x");
            u.setSelect({ left: 0, top: 0, width: 0, height: 0 }, false);
            setRange(a, b);
          }],
        },
      }, [[POS0, POS1], [null, null]], $("annplot"));
      annPlot.setScale("x", { min: x0, max: x1 });
      attachWheel(annPlot);
      attachPan(annPlot);
      wireAnnotation(annPlot);
    }
    attachWheel(main);
    attachPan(main);
    if (gcPlot) attachPan(gcPlot);

    refreshMain();

    $("reset").onclick = () => setRange(POS0, POS1);
    $("clip").onchange = refreshMain;
    $("fity").onchange = refreshMain;
    if (META.has_gc) {
      $("gcline").onchange = () => main.setSeries(ICOL_GC, { show: $("gcline").checked });
      $("gcpanel").onchange = () => { $("gcplot").style.display = $("gcpanel").checked ? "" : "none"; refreshMain(); };
      $("gcplot").style.display = "none"; // the GC curve is on the coverage plot (right axis); the panel is optional
    } else document.querySelectorAll(".gcopt").forEach((el) => (el.style.display = "none"));
    window.__viewer = { main, gcPlot, get annPlot() { return annPlot; }, ann: ANN }; // handles for the tests
    window.addEventListener("resize", () => {
      const w2 = widthOf(); main.setSize({ width: w2, height: 420 });
      if (annPlot) annPlot.setSize({ width: w2, height: META.ann.lanes * LANE_H + 8 });
      if (gcPlot) gcPlot.setSize({ width: w2, height: 130 });
      refreshMain();
    });
  }

  // ---- ROI table ----------------------------------------------------------
  const COLS = [["start", "Start"], ["end", "End"], ["size", "Size (bp)"], ["type", "Type"], ["mean_cov", LBL.mean_signal], ["mean_rm", LBL.mean_baseline], ["log2", "log2 ratio"], ["maxz", LBL.score]].concat(ANN ? [["geneStr", "Gene"], ["productStr", "Product"]] : []);
  let sortKey = "size", sortAsc = false;
  const esc = (t) => String(t).replace(/[&<>]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c]));
  function renderTable() {
    const minSize = +$("minsize").value || 0, type = $("type").value;
    const q = ANN ? $("search").value.trim().toLowerCase() : "";
    const rows = ROIS.filter((r) => r.size >= minSize && (type === "all" || r.type === type) &&
      (!q || r.genes.some((i) => `${ANN[i].name} ${ANN[i].product} ${ANN[i].id}`.toLowerCase().includes(q))));
    rows.sort((a, b) => {
      const x = a[sortKey], y = b[sortKey];
      const c = typeof x === "string" ? x.localeCompare(y) : (x - y) || 0;
      return sortAsc ? c : -c;
    });
    $("count").textContent = `${rows.length} of ${ROIS.length} regions`;
    const head = COLS.map(([k, label]) => `<th data-k="${k}" class="${(k === sortKey ? "sorted" + (sortAsc ? " asc" : "") : "") + (k === "geneStr" || k === "productStr" ? " txt" : "")}">${label}</th>`).join("");
    const num = (v, d) => (Number.isFinite(v) ? v.toFixed(d) : "-");
    const body = rows.map((r, i) => `<tr data-i="${ROIS.indexOf(r)}"><td>${fmt(r.start)}</td><td>${fmt(r.end)}</td><td>${fmt(r.size)}</td>` +
      `<td><span class="tag ${r.type}">${r.type === "low" ? LBL.low_short : LBL.high_short}</span></td><td>${num(r.mean_cov, 1)}</td><td>${num(r.mean_rm, 1)}</td><td>${num(r.log2, 2)}</td><td>${num(r.maxz, 1)}</td>` +
      (ANN ? `<td class="txt" title="${esc(r.geneStr)}">${esc(r.geneShort) || "<span style='color:#999'>intergenic</span>"}</td><td class="txt" title="${esc(r.productStr)}">${esc(r.productStr)}</td>` : "") + "</tr>").join("");
    $("roitable").innerHTML = `<thead><tr>${head}</tr></thead><tbody>${body}</tbody>`;
  }
  function wireTable() {
    $("roitable").addEventListener("click", (e) => {
      const th = e.target.closest("th");
      if (th) { const k = th.dataset.k; sortAsc = k === sortKey ? !sortAsc : false; sortKey = k; renderTable(); return; }
      const tr = e.target.closest("tr[data-i]");
      if (!tr) return;
      const r = ROIS[+tr.dataset.i];
      selected = r;
      const pad = Math.max(1000, r.size);
      setRange(r.start - pad, r.end + pad);
      main.redraw(false); if (annPlot) annPlot.redraw(false);
      document.querySelectorAll("tbody tr.sel").forEach((t) => t.classList.remove("sel"));
      tr.classList.add("sel");
      $("plot").scrollIntoView({ behavior: "smooth", block: "center" });
    });
    $("minsize").oninput = renderTable; $("type").onchange = renderTable;
    if (ANN) $("search").oninput = renderTable; else document.querySelectorAll(".annopt").forEach((el) => (el.style.display = "none"));
  }

  // ---- start --------------------------------------------------------------
  (async function () {
    try {
      status("decoding data...");
      LEVELS = await decodeLevels();
      build(); renderTable(); wireTable();
      window.__ready = true;
    } catch (err) { status("Error: " + err.message, true); throw err; }
  })();
})();
