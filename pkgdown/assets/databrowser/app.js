/* Data browser for the spacelink pkgdown site.
 *
 * Reads the static bundle written by data-raw/build_databrowser_data.R and draws
 * the spatial expression of a gene alongside cell-type abundance, plus the ESV
 * and PoPS readouts for the current selection. Everything runs client-side; the
 * only network traffic is the binary payload under ../databrowser/.
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

  var state = {
    manifest: null,
    datasetId: null,
    meta: null,
    coords: null,
    cellProp: null,
    esv: null,
    geneIdx: -1,
    cellType: "whole",
    diseaseIdx: -1,
    exprChunks: {},
    popsCache: {},
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
    '  <div class="dbx-field"><label for="dbx-celltype">Cell type</label>',
    '    <select id="dbx-celltype"></select></div>',
    '  <div class="dbx-field"><label for="dbx-disease">Disease (PoPS)</label>',
    '    <select id="dbx-disease"></select></div>',
    '</div>',
    '<div class="dbx-scores">',
    '  <div class="dbx-tile"><div class="dbx-tile-label">ESV score</div>',
    '    <div class="dbx-tile-value" id="dbx-esv">—</div>',
    '    <div class="dbx-tile-sub" id="dbx-esv-sub"></div></div>',
    '  <div class="dbx-tile"><div class="dbx-tile-label">PoPS score</div>',
    '    <div class="dbx-tile-value" id="dbx-pops">—</div>',
    '    <div class="dbx-tile-sub" id="dbx-pops-sub">No disease selected</div></div>',
    '</div>',
    '<div class="dbx-status" id="dbx-status">Loading…</div>',
    '<div class="dbx-plots dbx-two" id="dbx-plots">',
    '  <div class="dbx-panel" id="dbx-expr-panel">',
    '    <h3 id="dbx-expr-title">Gene expression</h3>',
    '    <div class="dbx-sub" id="dbx-expr-sub"></div>',
    '    <div class="dbx-canvas-wrap"><canvas id="dbx-expr-canvas"></canvas>',
    '      <div class="dbx-tip" id="dbx-expr-tip"></div></div>',
    '    <div class="dbx-legend"><span id="dbx-expr-lo">0</span>',
    '      <div class="dbx-bar" id="dbx-expr-bar"></div>',
    '      <span id="dbx-expr-hi"></span></div>',
    '    <details class="dbx-details"><summary>Values as a table</summary>',
    '      <div id="dbx-expr-table"></div></details></div>',
    '  <div class="dbx-panel" id="dbx-prop-panel">',
    '    <h3 id="dbx-prop-title">Cell type proportion</h3>',
    '    <div class="dbx-sub" id="dbx-prop-sub"></div>',
    '    <div id="dbx-prop-single"><div class="dbx-canvas-wrap">',
    '      <canvas id="dbx-prop-canvas"></canvas>',
    '      <div class="dbx-tip" id="dbx-prop-tip"></div></div></div>',
    '    <div id="dbx-prop-grid" class="dbx-mini-grid"></div>',
    '    <div class="dbx-legend"><span id="dbx-prop-lo">0</span>',
    '      <div class="dbx-bar" id="dbx-prop-bar"></div>',
    '      <span id="dbx-prop-hi"></span></div>',
    '    <details class="dbx-details"><summary>Values as a table</summary>',
    '      <div id="dbx-prop-table"></div></details></div>',
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
     ratio. The map is mirrored on both axes: the vertical flip matches the
     package vignettes (ggplot2 scale_y_reverse) and the horizontal flip puts
     the section in the requested orientation. */
  function paint(canvas, vals, vmax, rampName, cssW, maxH, minR) {
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

    var px = new Float32Array(n), py = new Float32Array(n);
    for (var j = 0; j < n; j++) {
      px[j] = ox + (xmax - coords[j * 2]) * s;        // x mirrored
      py[j] = oy + (ymax - coords[j * 2 + 1]) * s;    // y mirrored
    }

    // Radius from the typical spot spacing so dots read as a tissue section
    // rather than a sparse cloud, clamped to stay legible at any width.
    var r = Math.max(minR || 1.2,
      Math.min(6, 0.45 * Math.sqrt((dx * s) * (dy * s) / n)));

    var lut = buildLut(rampFor(rampName));
    var scale = vmax > 0 ? 255 / vmax : 0;
    for (var k = 0; k < n; k++) {
      var q = Math.max(0, Math.min(255, Math.round(vals[k] * scale)));
      ctx.fillStyle = "rgb(" + lut[q * 3] + "," + lut[q * 3 + 1] + "," + lut[q * 3 + 2] + ")";
      ctx.beginPath();
      ctx.arc(px[k], py[k], r, 0, 6.283185307179586);
      ctx.fill();
    }
    return { px: px, py: py, r: r };
  }

  function setLegend(prefix, vmax, rampName) {
    el("#" + prefix + "-bar").style.background =
      "linear-gradient(to right," + rampFor(rampName).join(",") + ")";
    el("#" + prefix + "-lo").textContent = "0";
    el("#" + prefix + "-hi").textContent = fmt(vmax, vmax < 1 ? 3 : 2);
  }

  function attachHover(canvas, tip, geom, vals, unitLabel, prefixLabel) {
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
        unitLabel + ": <b>" + fmt(vals[best], 3) + "</b>";
      tip.style.display = "block";
      var tw = tip.offsetWidth, th = tip.offsetHeight;
      tip.style.left = Math.max(0, Math.min(geom.px[best] + 12, canvas.clientWidth - tw)) + "px";
      tip.style.top = Math.max(0, geom.py[best] - th - 8) + "px";
    };
    canvas.onmouseleave = function () { tip.style.display = "none"; };
  }

  function summaryTable(target, vals, unitLabel) {
    var arr = Array.prototype.slice.call(vals).sort(function (a, b) { return a - b; });
    var sum = 0, nz = 0;
    for (var i = 0; i < vals.length; i++) { sum += vals[i]; if (vals[i] > 0) nz++; }
    var rows = [
      ["Minimum", fmt(arr[0], 3)],
      ["25th percentile", fmt(quantile(arr, 0.25), 3)],
      ["Median", fmt(quantile(arr, 0.5), 3)],
      ["75th percentile", fmt(quantile(arr, 0.75), 3)],
      ["Maximum", fmt(arr[arr.length - 1], 3)],
      ["Mean", fmt(sum / vals.length, 3)],
      ["Spots above zero", nz + " of " + vals.length]
    ];
    target.innerHTML = "<table><thead><tr><th scope=\"col\">Statistic</th><th scope=\"col\">" +
      esc(unitLabel) + "</th></tr></thead><tbody>" +
      rows.map(function (r) {
        return "<tr><th scope=\"row\">" + r[0] + "</th><td>" + r[1] + "</td></tr>";
      }).join("") + "</tbody></table>";
  }

  // ---- data access ------------------------------------------------------

  function geneValues(geneIdx) {
    var meta = state.meta;
    var chunk = Math.floor(geneIdx / meta.chunkSize);
    var key = state.datasetId + ":" + chunk;
    var p = state.exprChunks[key] || fetchBin(BASE + "/" + state.datasetId +
      "/expr_" + String(chunk).padStart(3, "0") + ".bin", Uint8Array);
    state.exprChunks[key] = p;
    return p.then(function (u8) {
      var off = (geneIdx - chunk * meta.chunkSize) * meta.nSpots;
      var vmax = meta.exprMax[geneIdx];
      var out = new Float32Array(meta.nSpots);
      for (var i = 0; i < meta.nSpots; i++) out[i] = (u8[off + i] / 255) * vmax;
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

  function popsValue(geneIdx, diseaseIdx) {
    var row = state.meta.popsRow[geneIdx];
    if (row === null || row === undefined || row < 0) return Promise.resolve(null);
    var key = String(diseaseIdx);
    var p = state.popsCache[key] ||
      fetchBin(BASE + "/pops/" + diseaseIdx + ".bin", Float32Array);
    state.popsCache[key] = p;
    return p.then(function (f32) { return f32[row]; });
  }

  // ---- rendering --------------------------------------------------------

  function renderScores() {
    var meta = state.meta;
    var gi = state.geneIdx;
    var ctIdx = meta.cellTypes.indexOf(state.cellType);
    var esv = state.esv[ctIdx * meta.nGenes + gi];
    var esvEl = el("#dbx-esv");
    var ok = isFinite(esv);
    esvEl.textContent = ok ? esv.toFixed(4) : "Not available";
    esvEl.className = "dbx-tile-value" + (ok ? "" : " dbx-na");
    el("#dbx-esv-sub").textContent = meta.genes[gi] + " · " + state.cellType;

    var popsEl = el("#dbx-pops");
    var popsSub = el("#dbx-pops-sub");
    if (state.diseaseIdx < 0) {
      popsEl.textContent = "—";
      popsEl.className = "dbx-tile-value dbx-na";
      popsSub.textContent = "No disease selected";
      return;
    }
    var disease = state.manifest.diseases[state.diseaseIdx];
    var gene = meta.genes[gi];
    var want = ++state.token;
    popsValue(gi, state.diseaseIdx).then(function (v) {
      if (want !== state.token) return;
      if (v === null) {
        popsEl.textContent = "Not in PoPS table";
        popsEl.className = "dbx-tile-value dbx-na";
        popsSub.textContent = gene + " has no PoPS score";
      } else {
        popsEl.textContent = v.toFixed(3);
        popsEl.className = "dbx-tile-value";
        popsSub.textContent = gene + " · " + disease;
      }
    });
  }

  function renderExpression() {
    var meta = state.meta;
    var gi = state.geneIdx;
    var gene = meta.genes[gi];
    el("#dbx-expr-title").textContent = gene + " expression";
    el("#dbx-expr-sub").textContent = "Normalised count per spot · " + meta.label;

    var canvas = el("#dbx-expr-canvas");
    var width = canvas.parentNode.clientWidth || 360;
    var want = state.token;
    return geneValues(gi).then(function (vals) {
      if (want !== state.token) return;
      var vmax = meta.exprMax[gi] || 1;
      var geom = paint(canvas, vals, vmax, "blue", width, 560);
      setLegend("dbx-expr", vmax, "blue");
      attachHover(canvas, el("#dbx-expr-tip"), geom, vals, "Normalised count");
      summaryTable(el("#dbx-expr-table"), vals, "Normalised count");
    });
  }

  /* One proportion map for the selected cell type. */
  function renderSingleProportion() {
    var meta = state.meta;
    el("#dbx-prop-single").style.display = "";
    el("#dbx-prop-grid").style.display = "none";

    var col = meta.cellTypes.indexOf(state.cellType) - 1;
    var vals = propValues(col);
    var pmax = 0;
    for (var i = 0; i < vals.length; i++) if (vals[i] > pmax) pmax = vals[i];

    el("#dbx-prop-title").textContent = state.cellType + " proportion";
    el("#dbx-prop-sub").textContent = "Estimated proportion per spot · " + meta.label;

    var canvas = el("#dbx-prop-canvas");
    var width = canvas.parentNode.clientWidth || 360;
    var geom = paint(canvas, vals, pmax || 1, "orange", width, 560);
    setLegend("dbx-prop", pmax || 1, "orange");
    attachHover(canvas, el("#dbx-prop-tip"), geom, vals, "Proportion");
    summaryTable(el("#dbx-prop-table"), vals, "Proportion");
  }

  /* Small multiples: every cell type at once, on one shared colour scale so the
     panels are comparable - a per-panel maximum would make a cell type that
     never exceeds 0.37 look as abundant as one that reaches 0.97. */
  function renderProportionGrid() {
    var meta = state.meta;
    var types = meta.cellTypes.slice(1);
    var grid = el("#dbx-prop-grid");
    el("#dbx-prop-single").style.display = "none";
    grid.style.display = "";

    var cols = types.map(function (_, i) { return propValues(i); });
    var shared = 0;
    cols.forEach(function (v) {
      for (var i = 0; i < v.length; i++) if (v[i] > shared) shared = v[i];
    });

    el("#dbx-prop-title").textContent = "Cell type proportions";
    el("#dbx-prop-sub").textContent =
      types.length + " cell types · shared colour scale · " + meta.label;

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

    // Table view: one row per cell type, so the grid is readable without colour.
    el("#dbx-prop-table").innerHTML =
      '<table><thead><tr><th scope="col">Cell type</th><th scope="col">Median</th>' +
      '<th scope="col">Maximum</th></tr></thead><tbody>' +
      types.map(function (t, i) {
        var arr = Array.prototype.slice.call(cols[i]).sort(function (a, b) { return a - b; });
        return "<tr><th scope=\"row\">" + esc(t) + "</th><td>" +
          fmt(quantile(arr, 0.5), 3) + "</td><td>" + fmt(arr[arr.length - 1], 3) + "</td></tr>";
      }).join("") + "</tbody></table>";
  }

  function renderPlots() {
    if (!state.meta || state.geneIdx < 0) return;
    state.token++;

    if (state.cellType === "whole") renderProportionGrid();
    else renderSingleProportion();

    renderExpression().then(function () { setStatus(""); }).catch(function (e) {
      setStatus("Could not load expression data: " + e.message, true);
    });

    renderScores();
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
    renderPlots();
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

  el("#dbx-celltype").addEventListener("change", function () {
    state.cellType = this.value;
    renderPlots();
  });

  el("#dbx-disease").addEventListener("change", function () {
    state.diseaseIdx = +this.value;
    renderScores();
  });

  el("#dbx-dataset").addEventListener("change", function () {
    loadDataset(this.value);
  });

  // Redraw on theme flip so the ramps re-anchor for the new surface.
  new MutationObserver(renderPlots)
    .observe(document.documentElement, { attributes: true, attributeFilter: ["data-bs-theme"] });

  var resizeTimer;
  window.addEventListener("resize", function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(renderPlots, 150);
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

      var ctSel = el("#dbx-celltype");
      ctSel.innerHTML = state.meta.cellTypes.map(function (c) {
        return '<option value="' + esc(c) + '">' + esc(c) + "</option>";
      }).join("");
      if (state.meta.cellTypes.indexOf(state.cellType) < 0) state.cellType = "whole";
      ctSel.value = state.cellType;

      el("#dbx-prop-grid").innerHTML = "";   // cell types may differ per dataset

      var genes = state.meta.genes;
      var preferred = ["SNAP25", "MBP", "PLP1", "GFAP"];
      var start = 0;
      for (var i = 0; i < preferred.length; i++) {
        var at = genes.indexOf(preferred[i]);
        if (at >= 0) { start = at; break; }
      }
      state.geneIdx = start;
      geneInput.value = genes[start];
      renderPlots();
    }).catch(function (e) {
      setStatus("Could not load dataset: " + e.message, true);
    });
  }

  Promise.all([fetchJson(BASE + "/manifest.json")]).then(function (res) {
    state.manifest = res[0];

    el("#dbx-dataset").innerHTML = state.manifest.datasets.map(function (d) {
      return '<option value="' + esc(d.id) + '">' + esc(d.label) + "</option>";
    }).join("");

    el("#dbx-disease").innerHTML = '<option value="-1">None</option>' +
      state.manifest.diseases.map(function (d, i) {
        return '<option value="' + i + '">' + esc(d) + "</option>";
      }).join("");

    return loadDataset(state.manifest.datasets[0].id);
  }).catch(function (e) {
    setStatus("Could not load the data browser: " + e.message, true);
  });
})();
