/* Data browser for the spacelink pkgdown site.
 *
 * Reads the static bundle written by data-raw/build_databrowser_data.R. Picking a
 * dataset and a gene draws that gene's spatial expression beside the abundance of
 * every cell type, and lists the gene's ESV scores and PoPS disease scores.
 * Everything runs client-side; the only network traffic is the binary payload
 * under ../databrowser/.
 */
(function () {
  "use strict";

  var BASE = "../databrowser";
  var mount = document.getElementById("spacelink-databrowser");
  if (!mount) return;

  // Sequential ramps: one hue each, light -> dark (project palette, blue steps
  // 100-700). Expression takes blue; the second sequential context on screen
  // takes orange. Dark mode flips the anchor so "near zero" recedes into the
  // dark surface instead of glowing against it.
  var RAMPS = {
    blue: ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"],
    orange: ["#ffd5c6", "#ffa98a", "#ff743e", "#e05104", "#b43f02", "#8a2e00", "#611e01"]
  };

  // Gene each dataset opens on: a marker whose spatial pattern reads clearly in
  // that tissue. A dataset missing from this map, or naming a gene its panel
  // does not carry, falls back to the first gene alphabetically.
  var DEFAULT_GENE = {
    visium_human_dlpfc: "SNAP25",
    cosmx_human_frontal_cortex: "SNAP25",
    visium_human_liver: "APOA1",
    cosmx_human_liver: "APOA1",
    visium_human_lymph_node: "IL7R",
    cosmx_human_lymph: "IL7R"
  };

  var state = {
    manifest: null,
    datasetId: null,
    meta: null,
    coords: null,
    cellProp: null,
    esv: null,
    pops: null,
    geneIdx: -1,
    exprChunks: {},
    token: 0
  };

  // ---- small helpers ----------------------------------------------------

  function el(sel) { return mount.querySelector(sel); }

  function hexToRgb(h) {
    return [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
  }

  /* 256-entry lookup built by interpolating between ramp stops. */
  function buildLut(stops) {
    var rgb = stops.map(hexToRgb);
    var lut = new Uint8Array(256 * 3);
    for (var i = 0; i < 256; i++) {
      var t = (i / 255) * (rgb.length - 1);
      var k = Math.min(Math.floor(t), rgb.length - 2);
      var f = t - k;
      for (var c = 0; c < 3; c++) {
        lut[i * 3 + c] = Math.round(rgb[k][c] + (rgb[k + 1][c] - rgb[k][c]) * f);
      }
    }
    return lut;
  }

  function isDark() {
    return document.documentElement.getAttribute("data-bs-theme") === "dark";
  }

  function rampFor(name) {
    var stops = RAMPS[name];
    return isDark() ? stops.slice().reverse() : stops;
  }

  function fmt(v, d) {
    if (v === null || v === undefined || !isFinite(v)) return "—";
    return v.toFixed(d === undefined ? 3 : d);
  }

  /* Counts are integers, but quantiles and means of them are not; show a whole
     number when the value is one and a single decimal otherwise. */
  function fmtCount(v) {
    if (v === null || v === undefined || !isFinite(v)) return "—";
    return Number.isInteger(v) ? String(v) : v.toFixed(1);
  }

  function fmtInt(v) {
    return String(v).replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  }

  /* Colour-bar ticks span exact counts in the hundreds and binned means below
     one, so pick the precision from the magnitude rather than fixing it. */
  function fmtScale(v) {
    if (v === null || v === undefined || !isFinite(v)) return "—";
    if (Number.isInteger(v)) return String(v);
    if (v >= 10) return v.toFixed(1);
    if (v >= 1) return v.toFixed(2);
    return v.toFixed(3);
  }

  function fetchBin(url, ctor) {
    return fetch(url).then(function (r) {
      if (!r.ok) throw new Error(url + " (" + r.status + ")");
      return r.arrayBuffer();
    }).then(function (b) { return new ctor(b); });
  }

  function fetchJson(url) {
    return fetch(url).then(function (r) {
      if (!r.ok) throw new Error(url + " (" + r.status + ")");
      return r.json();
    });
  }

  function quantile(sorted, p) {
    if (!sorted.length) return NaN;
    var i = (sorted.length - 1) * p;
    var lo = Math.floor(i), hi = Math.ceil(i);
    return sorted[lo] + (sorted[hi] - sorted[lo]) * (i - lo);
  }

  function esc(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

  // ---- markup -----------------------------------------------------------

  mount.innerHTML = [
    '<div class="dbx-controls">',
    '  <div class="dbx-field"><label for="dbx-dataset">Dataset</label>',
    '    <select id="dbx-dataset"></select></div>',
    '  <div class="dbx-field dbx-combo"><label for="dbx-gene">Gene</label>',
    '    <input id="dbx-gene" type="text" autocomplete="off" spellcheck="false"',
    '      role="combobox" aria-autocomplete="list" aria-expanded="false"',
    '      aria-controls="dbx-gene-options" placeholder="Type to search…">',
    '    <div class="dbx-options" id="dbx-gene-options" role="listbox"></div></div>',
    '</div>',
    '<div class="dbx-status" id="dbx-status">Loading…</div>',
    '<div class="dbx-plots dbx-two" id="dbx-plots">',
    '  <div class="dbx-panel" id="dbx-expr-panel">',
    '    <h3 id="dbx-expr-title">Gene expression</h3>',
    '    <div class="dbx-sub" id="dbx-expr-sub"></div>',
    '    <div class="dbx-canvas-wrap"><canvas id="dbx-expr-canvas"></canvas>',
    '      <div class="dbx-tip" id="dbx-expr-tip"></div></div>',
    '    <div class="dbx-legend"><div class="dbx-bar" id="dbx-expr-bar"></div></div>',
    '    <div class="dbx-ticks"><span id="dbx-expr-lo">0</span>',
    '      <span id="dbx-expr-mid"></span><span id="dbx-expr-hi"></span></div>',
    '    <div class="dbx-legend-note" id="dbx-expr-note"></div>',
    '    <details class="dbx-details"><summary>Values as a table</summary>',
    '      <div id="dbx-expr-table"></div></details></div>',
    '  <div class="dbx-panel" id="dbx-prop-panel">',
    '    <h3 id="dbx-prop-title">Cell type proportions</h3>',
    '    <div class="dbx-sub" id="dbx-prop-sub"></div>',
    '    <div id="dbx-prop-grid" class="dbx-mini-grid"></div>',
    '    <div class="dbx-legend"><div class="dbx-bar" id="dbx-prop-bar"></div></div>',
    '    <div class="dbx-ticks"><span id="dbx-prop-lo">0</span>',
    '      <span id="dbx-prop-mid"></span><span id="dbx-prop-hi"></span></div>',
    '    <details class="dbx-details"><summary>Values as a table</summary>',
    '      <div id="dbx-prop-table"></div></details></div>',
    '</div>',
    '<div class="dbx-tables">',
    '  <div class="dbx-panel"><h3>ESV scores</h3>',
    '    <div class="dbx-sub" id="dbx-esv-sub"></div>',
    '    <div class="dbx-scroll" id="dbx-esv-table"></div></div>',
    '  <div class="dbx-panel"><h3>PoPS disease scores</h3>',
    '    <div class="dbx-sub" id="dbx-pops-sub"></div>',
    '    <div class="dbx-scroll" id="dbx-pops-table"></div></div>',
    '</div>'
  ].join("\n");

  var status = el("#dbx-status");

  function setStatus(msg, isError) {
    status.textContent = msg || "";
    status.style.display = msg ? "" : "none";
    status.className = "dbx-status" + (isError ? " dbx-error" : "");
  }

  // ---- geometry & painting ----------------------------------------------

  /* Projects the spot coordinates into a cssW-wide canvas, preserving aspect
     ratio. Both axes run in the same direction as the stored coordinates.
     `logScale` maps colour by log1p, which raw counts need: their distribution
     is steep enough that a linear ramp puts ~85% of spots in the lowest band
     and the tissue structure disappears. */
  function paint(canvas, vals, vmax, rampName, cssW, maxH, minR, logScale) {
    var coords = state.coords;
    var n = state.meta.nSpots;

    var xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (var i = 0; i < n; i++) {
      var x = coords[i * 2], y = coords[i * 2 + 1];
      if (x < xmin) xmin = x; if (x > xmax) xmax = x;
      if (y < ymin) ymin = y; if (y > ymax) ymax = y;
    }
    var dx = (xmax - xmin) || 1, dy = (ymax - ymin) || 1;
    var pad = Math.max(3, Math.round(cssW * 0.025));

    var cssH = Math.round((cssW - 2 * pad) * (dy / dx)) + 2 * pad;
    cssH = Math.max(60, Math.min(cssH, maxH));

    var dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
    canvas.style.width = cssW + "px";
    canvas.style.height = cssH + "px";

    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);

    var s = Math.min((cssW - 2 * pad) / dx, (cssH - 2 * pad) / dy);
    var ox = (cssW - dx * s) / 2, oy = (cssH - dy * s) / 2;

    // CosMx stores its y axis in the opposite direction to Visium, so those
    // sections are mirrored vertically to match the orientation used elsewhere.
    var flipY = state.meta.flipY !== undefined
      ? state.meta.flipY
      : /^cosmx_/.test(state.datasetId);   // fallback for bundles predating flipY

    var px = new Float32Array(n), py = new Float32Array(n);
    for (var j = 0; j < n; j++) {
      px[j] = ox + (coords[j * 2] - xmin) * s;
      py[j] = flipY ? oy + (ymax - coords[j * 2 + 1]) * s
                    : oy + (coords[j * 2 + 1] - ymin) * s;
    }

    // Radius from the typical spot spacing so dots read as a tissue section
    // rather than a sparse cloud, clamped to stay legible at any width.
    var r = Math.max(minR || 1.2,
      Math.min(6, 0.45 * Math.sqrt((dx * s) * (dy * s) / n)));

    var lut = buildLut(rampFor(rampName));
    var denom = logScale ? Math.log1p(vmax) : vmax;
    var scale = denom > 0 ? 255 / denom : 0;

    // Bucket the points by colour first. A binned CosMx panel draws ~20k marks
    // across 14 canvases, and changing fillStyle per mark dominates the cost;
    // this caps it at 256 style changes and one path per colour.
    var buckets = new Array(256);
    for (var k = 0; k < n; k++) {
      var t = logScale ? Math.log1p(vals[k]) : vals[k];
      var q = Math.max(0, Math.min(255, Math.round(t * scale)));
      (buckets[q] || (buckets[q] = [])).push(k);
    }
    var square = r < 2;          // below ~2px a rect and a disc are identical
    var d = r * 2;
    for (var qi = 0; qi < 256; qi++) {
      var b = buckets[qi];
      if (!b) continue;
      ctx.fillStyle = "rgb(" + lut[qi * 3] + "," + lut[qi * 3 + 1] + "," + lut[qi * 3 + 2] + ")";
      if (square) {
        for (var j = 0; j < b.length; j++) ctx.fillRect(px[b[j]] - r, py[b[j]] - r, d, d);
      } else {
        ctx.beginPath();
        for (var j2 = 0; j2 < b.length; j2++) {
          ctx.moveTo(px[b[j2]] + r, py[b[j2]]);
          ctx.arc(px[b[j2]], py[b[j2]], r, 0, 6.283185307179586);
        }
        ctx.fill();
      }
    }
    return { px: px, py: py, r: r };
  }

  /* On a log ramp the bar's midpoint is not vmax/2, so label it explicitly -
     without that tick the gradient reads as linear and overstates the low end. */
  function setLegend(prefix, vmax, rampName, logScale) {
    el("#" + prefix + "-bar").style.background =
      "linear-gradient(to right," + rampFor(rampName).join(",") + ")";
    el("#" + prefix + "-lo").textContent = "0";
    el("#" + prefix + "-hi").textContent = logScale ? fmtScale(vmax) : fmt(vmax, vmax < 1 ? 3 : 2);
    var mid = el("#" + prefix + "-mid");
    if (mid) {
      if (logScale) {
        var midVal = Math.expm1(0.5 * Math.log1p(vmax));
        // Whole counts read better rounded; sub-unit bin means need decimals.
        mid.textContent = vmax >= 10 ? String(Math.round(midVal)) : fmtScale(midVal);
        mid.style.display = "";
      } else {
        mid.style.display = "none";
      }
    }
  }

  function attachHover(canvas, tip, geom, vals, unitLabel, prefixLabel, isCount) {
    var hitR = Math.max(geom.r + 3, 6);
    canvas.onmousemove = function (ev) {
      var rect = canvas.getBoundingClientRect();
      var mx = ev.clientX - rect.left, my = ev.clientY - rect.top;
      var best = -1, bestD = hitR * hitR;
      for (var i = 0; i < geom.px.length; i++) {
        var ddx = geom.px[i] - mx, ddy = geom.py[i] - my;
        var d = ddx * ddx + ddy * ddy;
        if (d < bestD) { bestD = d; best = i; }
      }
      if (best < 0) { tip.style.display = "none"; return; }
      tip.innerHTML = (prefixLabel ? esc(prefixLabel) + "<br>" : "") +
        unitLabel + ": <b>" + (isCount ? String(vals[best]) : fmt(vals[best], 3)) + "</b>";
      tip.style.display = "block";
      var tw = tip.offsetWidth, th = tip.offsetHeight;
      tip.style.left = Math.max(0, Math.min(geom.px[best] + 12, canvas.clientWidth - tw)) + "px";
      tip.style.top = Math.max(0, geom.py[best] - th - 8) + "px";
    };
    canvas.onmouseleave = function () { tip.style.display = "none"; };
  }

  function summaryTable(target, vals, unitLabel, isCount) {
    var arr = Array.prototype.slice.call(vals).sort(function (a, b) { return a - b; });
    var sum = 0, nz = 0;
    for (var i = 0; i < vals.length; i++) { sum += vals[i]; if (vals[i] > 0) nz++; }
    var f = isCount ? fmtCount : function (v) { return fmt(v, 3); };
    var rows = [
      ["Minimum", f(arr[0])],
      ["25th percentile", f(quantile(arr, 0.25))],
      ["Median", f(quantile(arr, 0.5))],
      ["75th percentile", f(quantile(arr, 0.75))],
      ["Maximum", f(arr[arr.length - 1])],
      ["Mean", isCount ? (sum / vals.length).toFixed(2) : fmt(sum / vals.length, 3)],
      ["Spots above zero", nz + " of " + vals.length]
    ];
    target.innerHTML = "<table><thead><tr><th scope=\"col\">Statistic</th><th scope=\"col\">" +
      esc(unitLabel) + "</th></tr></thead><tbody>" +
      rows.map(function (r) {
        return "<tr><th scope=\"row\">" + r[0] + "</th><td>" + r[1] + "</td></tr>";
      }).join("") + "</tbody></table>";
  }

  /* Score table with an inline magnitude bar. The bar is a neutral tint, not a
     ramp colour - length carries the value, and no hue is spent implying a link
     to the maps above. */
  function scoreTable(target, rowsIn, keyHeader, valHeader, decimals) {
    var vmax = 0;
    rowsIn.forEach(function (r) { if (isFinite(r[1]) && r[1] > vmax) vmax = r[1]; });
    target.innerHTML = '<table class="dbx-score-table"><thead><tr><th scope="col">' +
      esc(keyHeader) + '</th><th scope="col">' + esc(valHeader) +
      "</th></tr></thead><tbody>" +
      rowsIn.map(function (r) {
        var ok = isFinite(r[1]);
        var pct = ok && vmax > 0 ? (r[1] / vmax) * 100 : 0;
        return '<tr><th scope="row">' + esc(r[0]) + "</th>" +
          '<td><span class="dbx-barcell" style="--w:' + pct.toFixed(1) + '%"></span>' +
          "<span>" + (ok ? r[1].toFixed(decimals) : "—") + "</span></td></tr>";
      }).join("") + "</tbody></table>";
  }

  // ---- data access ------------------------------------------------------

  function geneValues(geneIdx) {
    var meta = state.meta;
    var chunk = Math.floor(geneIdx / meta.chunkSize);
    var key = state.datasetId + ":" + chunk;
    var p = state.exprChunks[key] || fetchBin(BASE + "/" + state.datasetId +
      "/expr_" + String(chunk).padStart(3, "0") + ".bin", Uint16Array);
    state.exprChunks[key] = p;
    return p.then(function (u16) {
      var off = (geneIdx - chunk * meta.chunkSize) * meta.nSpots;
      var scale = meta.exprScale ? meta.exprScale[geneIdx] : 1;
      // Per-spot datasets store exact integers (scale 1); binned datasets store
      // fractional bin means scaled to each gene's own maximum.
      if (scale === 1) return u16.subarray(off, off + meta.nSpots);
      var out = new Float32Array(meta.nSpots);
      for (var i = 0; i < meta.nSpots; i++) out[i] = u16[off + i] * scale;
      return out;
    });
  }

  /* Column `col` of the spot-major cell-type proportion matrix. */
  function propValues(col) {
    var meta = state.meta;
    var nCT = meta.cellTypes.length - 1;      // cellProp has no "whole" column
    var out = new Float32Array(meta.nSpots);
    for (var i = 0; i < meta.nSpots; i++) out[i] = state.cellProp[i * nCT + col];
    return out;
  }

  // ---- rendering --------------------------------------------------------

  function renderExpression() {
    var meta = state.meta;
    var gi = state.geneIdx;
    var unit = meta.binned ? "Mean count" : "Count";
    el("#dbx-expr-title").textContent = meta.genes[gi] + " expression";
    el("#dbx-expr-sub").textContent = meta.binned
      ? "Mean raw count per binned cell · " + meta.label
      : "Raw count per spot · " + meta.label;

    var canvas = el("#dbx-expr-canvas");
    var width = canvas.parentNode.clientWidth || 360;
    var want = state.token;
    return geneValues(gi).then(function (vals) {
      if (want !== state.token) return;
      var vmax = meta.exprMax[gi] || 1;
      var geom = paint(canvas, vals, vmax, "blue", width, 560, undefined, true);
      setLegend("dbx-expr", vmax, "blue", true);
      el("#dbx-expr-note").textContent = meta.binned
        ? "Colour on a log scale. " + fmtInt(meta.nSource) + " cells aggregated into " +
          fmtInt(meta.nSpots) + " spatial bins; each bin shows the mean of its cells."
        : "Colour on a log scale; counts are exact.";
      attachHover(canvas, el("#dbx-expr-tip"), geom, vals, unit, null, !meta.binned);
      summaryTable(el("#dbx-expr-table"), vals, unit, !meta.binned);
    });
  }

  /* Small multiples: every cell type at once. They share one colour scale so the
     panels stay comparable - a per-panel maximum would make a cell type that
     never exceeds 0.37 look as abundant as one that reaches 0.97. */
  function renderProportionGrid() {
    var meta = state.meta;
    var grid = el("#dbx-prop-grid");

    // Fixed dataset order, deliberately not the ESV ranking used by the table:
    // a tile keeps its position as the gene changes, so a cell type can be
    // followed across genes.
    var types = meta.cellTypes.slice(1);
    var cols = types.map(function (_, i) { return propValues(i); });
    var shared = 0;
    cols.forEach(function (v) {
      for (var i = 0; i < v.length; i++) if (v[i] > shared) shared = v[i];
    });

    el("#dbx-prop-sub").textContent = types.length + " cell types · " +
      (meta.binned ? "mean over each bin · " : "") + meta.label;

    if (grid.childElementCount !== types.length) {
      grid.innerHTML = types.map(function (t, i) {
        return '<figure class="dbx-mini"><figcaption title="' + esc(t) + '">' +
          esc(t) + '</figcaption><div class="dbx-canvas-wrap">' +
          '<canvas id="dbx-mini-' + i + '"></canvas>' +
          '<div class="dbx-tip" id="dbx-mini-tip-' + i + '"></div></div></figure>';
      }).join("");
    }

    var cellW = grid.firstElementChild
      ? grid.firstElementChild.querySelector(".dbx-canvas-wrap").clientWidth
      : 120;
    if (!cellW) cellW = 120;

    types.forEach(function (t, i) {
      var canvas = el("#dbx-mini-" + i);
      var geom = paint(canvas, cols[i], shared || 1, "orange", cellW, 260, 0.8);
      attachHover(canvas, el("#dbx-mini-tip-" + i), geom, cols[i], "Proportion", t);
    });

    setLegend("dbx-prop", shared || 1, "orange");

    el("#dbx-prop-table").innerHTML =
      '<table><thead><tr><th scope="col">Cell type</th><th scope="col">Median</th>' +
      '<th scope="col">Maximum</th></tr></thead><tbody>' +
      types.map(function (t, i) {
        var arr = Array.prototype.slice.call(cols[i]).sort(function (a, b) { return a - b; });
        return "<tr><th scope=\"row\">" + esc(t) + "</th><td>" +
          fmt(quantile(arr, 0.5), 3) + "</td><td>" + fmt(arr[arr.length - 1], 3) + "</td></tr>";
      }).join("") + "</tbody></table>";
  }

  /* Every ESV score for this gene: the whole-tissue score first, then each
     cell type, in the order the dataset stores them. */
  /* Cell types ranked by the selected gene's ESV, strongest first, as indices
     into meta.cellTypes. Index 0 ("whole") is excluded - it is the whole-tissue
     score, not a cell type. Pairs without a score sort last, keeping their
     dataset order among themselves. Used by the ESV table only; the proportion
     maps stay in dataset order so tiles do not move between genes. */
  function cellTypeOrder() {
    var meta = state.meta, gi = state.geneIdx;
    var idx = [];
    for (var i = 1; i < meta.cellTypes.length; i++) idx.push(i);
    idx.sort(function (a, b) {
      var va = state.esv[a * meta.nGenes + gi];
      var vb = state.esv[b * meta.nGenes + gi];
      var na = !isFinite(va), nb = !isFinite(vb);
      if (na && nb) return a - b;
      if (na) return 1;
      if (nb) return -1;
      if (vb !== va) return vb - va;
      return a - b;
    });
    return idx;
  }

  function renderEsvTable() {
    var meta = state.meta;
    var gi = state.geneIdx;
    el("#dbx-esv-sub").textContent = meta.genes[gi] + " · " + meta.label;
    // "whole" stays pinned at the top; the cell types follow in ESV order.
    var rows = [["whole", state.esv[gi]]];
    cellTypeOrder().forEach(function (i) {
      rows.push([meta.cellTypes[i], state.esv[i * meta.nGenes + gi]]);
    });
    scoreTable(el("#dbx-esv-table"), rows, "Cell type", "ESV", 4);
  }

  /* Every PoPS score for this gene, strongest first - the useful read is which
     traits the gene is prioritised for. */
  function renderPopsTable() {
    var meta = state.meta;
    var gi = state.geneIdx;
    var gene = meta.genes[gi];
    var target = el("#dbx-pops-table");
    var row = meta.popsRow[gi];

    if (row === null || row === undefined || row < 0 || !state.pops) {
      el("#dbx-pops-sub").textContent = gene + " is not in the PoPS tables";
      target.innerHTML = '<p class="dbx-empty-note">No PoPS scores are available for ' +
        esc(gene) + ".</p>";
      return;
    }
    var diseases = state.manifest.diseases;
    var nD = diseases.length;
    var rows = diseases.map(function (d, k) { return [d, state.pops[row * nD + k]]; });
    rows.sort(function (a, b) { return b[1] - a[1]; });
    el("#dbx-pops-sub").textContent =
      gene + " · " + nD + " traits, highest first";
    scoreTable(target, rows, "Trait", "PoPS", 3);
  }

  function renderAll() {
    if (!state.meta || state.geneIdx < 0) return;
    state.token++;
    renderProportionGrid();
    renderEsvTable();
    renderPopsTable();
    renderExpression().then(function () { setStatus(""); }).catch(function (e) {
      setStatus("Could not load expression data: " + e.message, true);
    });
  }

  // ---- gene combobox ----------------------------------------------------

  var geneInput = el("#dbx-gene");
  var geneOptions = el("#dbx-gene-options");
  var activeIdx = -1;
  var shown = [];

  function closeOptions() {
    geneOptions.classList.remove("open");
    geneInput.setAttribute("aria-expanded", "false");
    activeIdx = -1;
  }

  function openOptions(list) {
    shown = list;
    if (!list.length) {
      geneOptions.innerHTML = '<div class="dbx-empty">No matching gene</div>';
    } else {
      geneOptions.innerHTML = list.map(function (g, i) {
        return '<div role="option" data-i="' + i + '" aria-selected="false">' + esc(g) + "</div>";
      }).join("");
    }
    geneOptions.classList.add("open");
    geneInput.setAttribute("aria-expanded", "true");
    activeIdx = -1;
  }

  function highlight(i) {
    var kids = geneOptions.children;
    if (activeIdx >= 0 && kids[activeIdx]) kids[activeIdx].setAttribute("aria-selected", "false");
    activeIdx = i;
    if (i >= 0 && kids[i]) {
      kids[i].setAttribute("aria-selected", "true");
      kids[i].scrollIntoView({ block: "nearest" });
    }
  }

  /* Prefix matches first - typing "GFA" should surface GFAP before a substring
     hit elsewhere - then substring matches, capped so the list stays usable. */
  function searchGenes(q) {
    var genes = state.meta.genes;
    if (!q) return genes.slice(0, 100);
    q = q.toUpperCase();
    var pre = [], sub = [];
    for (var i = 0; i < genes.length; i++) {
      var u = genes[i].toUpperCase();
      var at = u.indexOf(q);
      if (at === 0) pre.push(genes[i]);
      else if (at > 0) sub.push(genes[i]);
      if (pre.length >= 100) break;
    }
    return pre.concat(sub).slice(0, 100);
  }

  function chooseGene(name) {
    var i = state.meta.genes.indexOf(name);
    if (i < 0) return false;
    state.geneIdx = i;
    geneInput.value = name;
    closeOptions();
    renderAll();
    return true;
  }

  geneInput.addEventListener("input", function () { openOptions(searchGenes(geneInput.value.trim())); });
  geneInput.addEventListener("focus", function () { openOptions(searchGenes(geneInput.value.trim())); });

  geneInput.addEventListener("keydown", function (ev) {
    if (ev.key === "ArrowDown" || ev.key === "ArrowUp") {
      ev.preventDefault();
      if (!geneOptions.classList.contains("open")) openOptions(searchGenes(geneInput.value.trim()));
      if (!shown.length) return;
      highlight(ev.key === "ArrowDown"
        ? (activeIdx + 1) % shown.length
        : (activeIdx <= 0 ? shown.length - 1 : activeIdx - 1));
    } else if (ev.key === "Enter") {
      ev.preventDefault();
      chooseGene(activeIdx >= 0 ? shown[activeIdx] : (shown[0] || geneInput.value.trim()));
    } else if (ev.key === "Escape") {
      closeOptions();
    }
  });

  geneOptions.addEventListener("mousedown", function (ev) {
    var t = ev.target.closest("[data-i]");
    if (!t) return;
    ev.preventDefault();
    chooseGene(shown[+t.dataset.i]);
  });

  geneInput.addEventListener("blur", function () {
    setTimeout(function () {
      closeOptions();
      // Snap back to the last valid gene if the box was left half-typed.
      if (state.geneIdx >= 0) geneInput.value = state.meta.genes[state.geneIdx];
    }, 120);
  });

  // ---- selection wiring -------------------------------------------------

  el("#dbx-dataset").addEventListener("change", function () { loadDataset(this.value); });

  // Redraw on theme flip so the ramps re-anchor for the new surface.
  new MutationObserver(renderAll)
    .observe(document.documentElement, { attributes: true, attributeFilter: ["data-bs-theme"] });

  var resizeTimer;
  window.addEventListener("resize", function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(renderAll, 150);
  });

  // ---- boot -------------------------------------------------------------

  function loadDataset(id) {
    setStatus("Loading dataset…");
    state.datasetId = id;
    var dir = BASE + "/" + id;
    return Promise.all([
      fetchJson(dir + "/meta.json"),
      fetchBin(dir + "/coords.bin", Float32Array),
      fetchBin(dir + "/celltype.bin", Float32Array),
      fetchBin(dir + "/esv.bin", Float32Array)
    ]).then(function (res) {
      state.meta = res[0];
      state.coords = res[1];
      state.cellProp = res[2];
      state.esv = res[3];
      state.exprChunks = {};
      el("#dbx-prop-grid").innerHTML = "";   // cell types may differ per dataset

      var genes = state.meta.genes;
      var start = genes.indexOf(DEFAULT_GENE[id] || "");
      if (start < 0) start = 0;      // unknown dataset, or gene absent from it
      state.geneIdx = start;
      geneInput.value = genes[start];
      renderAll();
    }).catch(function (e) {
      setStatus("Could not load dataset: " + e.message, true);
    });
  }

  fetchJson(BASE + "/manifest.json").then(function (manifest) {
    state.manifest = manifest;
    el("#dbx-dataset").innerHTML = manifest.datasets.map(function (d) {
      return '<option value="' + esc(d.id) + '">' + esc(d.label) + "</option>";
    }).join("");

    // PoPS is shared across datasets and small enough to fetch once up front.
    fetchBin(BASE + "/pops.bin", Float32Array).then(function (p) {
      state.pops = p;
      if (state.meta) renderPopsTable();
    }).catch(function () { state.pops = null; });

    return loadDataset(manifest.datasets[0].id);
  }).catch(function (e) {
    setStatus("Could not load the data browser: " + e.message, true);
  });
})();
