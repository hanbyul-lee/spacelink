/* Data browser for the spacelink pkgdown site.
 *
 * Reads the static bundle written by data-raw/build_databrowser_data.R and draws
 * two spatial scatter plots (gene expression, cell-type proportion) plus the ESV
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
    popsGenes: null,
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
    '<div class="dbx-plots" id="dbx-plots"></div>'
  ].join("\n");

  var status = el("#dbx-status");

  function setStatus(msg, isError) {
    status.textContent = msg || "";
    status.style.display = msg ? "" : "none";
    status.className = "dbx-status" + (isError ? " dbx-error" : "");
  }

  // ---- panels -----------------------------------------------------------

  function makePanel(id, title) {
    var d = document.createElement("div");
    d.className = "dbx-panel";
    d.innerHTML = [
      '<h3 id="' + id + '-title">' + title + "</h3>",
      '<div class="dbx-sub" id="' + id + '-sub"></div>',
      '<div class="dbx-canvas-wrap"><canvas id="' + id + '-canvas"></canvas>',
      '  <div class="dbx-tip" id="' + id + '-tip"></div></div>',
      '<div class="dbx-legend"><span id="' + id + '-lo">0</span>',
      '  <div class="dbx-bar" id="' + id + '-bar"></div>',
      '  <span id="' + id + '-hi"></span></div>',
      '<details class="dbx-details"><summary>Values as a table</summary>',
      '  <div id="' + id + '-table"></div></details>'
    ].join("\n");
    return d;
  }

  var plots = el("#dbx-plots");
  plots.appendChild(makePanel("dbx-expr", "Gene expression"));
  plots.appendChild(makePanel("dbx-prop", "Cell type proportion"));

  // ---- drawing ----------------------------------------------------------

  /* Draws one spatial scatter. `vals` is per-spot; `vmax` anchors the ramp at
     its top so the colour bar labels match. */
  function draw(panelId, vals, vmax, rampName, unitLabel) {
    var canvas = el("#" + panelId + "-canvas");
    var wrap = canvas.parentNode;
    var cssW = wrap.clientWidth || 360;
    var coords = state.coords;
    var n = state.meta.nSpots;

    // Data extent, y reversed to match the tissue orientation used in the
    // package vignettes (ggplot2 scale_y_reverse).
    var xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (var i = 0; i < n; i++) {
      var x = coords[i * 2], y = coords[i * 2 + 1];
      if (x < xmin) xmin = x; if (x > xmax) xmax = x;
      if (y < ymin) ymin = y; if (y > ymax) ymax = y;
    }
    var dx = (xmax - xmin) || 1, dy = (ymax - ymin) || 1;
    var pad = 10;
    var cssH = Math.round((cssW - 2 * pad) * (dy / dx)) + 2 * pad;
    cssH = Math.max(160, Math.min(cssH, 560));

    var dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
    canvas.style.height = cssH + "px";

    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);

    // Equal aspect: one scale for both axes, centred.
    var s = Math.min((cssW - 2 * pad) / dx, (cssH - 2 * pad) / dy);
    var ox = (cssW - dx * s) / 2, oy = (cssH - dy * s) / 2;

    var px = new Float32Array(n), py = new Float32Array(n);
    for (var j = 0; j < n; j++) {
      px[j] = ox + (coords[j * 2] - xmin) * s;
      py[j] = oy + (ymax - coords[j * 2 + 1]) * s;   // y reversed
    }

    // Radius from the typical spot spacing so dots read as a tissue section
    // rather than a sparse cloud, clamped to stay legible at any width.
    var r = Math.max(1.2, Math.min(6, 0.45 * Math.sqrt((dx * s) * (dy * s) / n)));

    var lut = buildLut(rampFor(rampName));
    var scale = vmax > 0 ? 255 / vmax : 0;
    for (var k = 0; k < n; k++) {
      var q = Math.max(0, Math.min(255, Math.round(vals[k] * scale)));
      ctx.fillStyle = "rgb(" + lut[q * 3] + "," + lut[q * 3 + 1] + "," + lut[q * 3 + 2] + ")";
      ctx.beginPath();
      ctx.arc(px[k], py[k], r, 0, 6.283185307179586);
      ctx.fill();
    }

    // colour bar + range labels
    el("#" + panelId + "-bar").style.background =
      "linear-gradient(to right," + rampFor(rampName).join(",") + ")";
    el("#" + panelId + "-lo").textContent = "0";
    el("#" + panelId + "-hi").textContent = fmt(vmax, vmax < 1 ? 3 : 2);

    renderTable(panelId, vals, unitLabel);
    attachHover(panelId, px, py, vals, r, unitLabel);
  }

  function renderTable(panelId, vals, unitLabel) {
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
    el("#" + panelId + "-table").innerHTML =
      '<table><caption class="visually-hidden">' + unitLabel +
      ' summary across spots</caption><thead><tr><th scope="col">Statistic</th>' +
      '<th scope="col">' + unitLabel + "</th></tr></thead><tbody>" +
      rows.map(function (r) {
        return "<tr><th scope=\"row\">" + r[0] + "</th><td>" + r[1] + "</td></tr>";
      }).join("") + "</tbody></table>";
  }

  function attachHover(panelId, px, py, vals, r, unitLabel) {
    var canvas = el("#" + panelId + "-canvas");
    var tip = el("#" + panelId + "-tip");
    var hitR = Math.max(r + 3, 6);

    canvas.onmousemove = function (ev) {
      var rect = canvas.getBoundingClientRect();
      var mx = ev.clientX - rect.left, my = ev.clientY - rect.top;
      var best = -1, bestD = hitR * hitR;
      for (var i = 0; i < px.length; i++) {
        var ddx = px[i] - mx, ddy = py[i] - my;
        var d = ddx * ddx + ddy * ddy;
        if (d < bestD) { bestD = d; best = i; }
      }
      if (best < 0) { tip.style.display = "none"; return; }
      tip.innerHTML = unitLabel + ": <b>" + fmt(vals[best], 3) + "</b>";
      tip.style.display = "block";
      var tw = tip.offsetWidth, th = tip.offsetHeight;
      tip.style.left = Math.max(0, Math.min(px[best] + 12, canvas.clientWidth - tw)) + "px";
      tip.style.top = Math.max(0, py[best] - th - 8) + "px";
    };
    canvas.onmouseleave = function () { tip.style.display = "none"; };
  }

  // ---- data access ------------------------------------------------------

  function geneValues(geneIdx) {
    var meta = state.meta;
    var chunk = Math.floor(geneIdx / meta.chunkSize);
    var key = state.datasetId + ":" + chunk;
    var cached = state.exprChunks[key];
    var p = cached || fetchBin(BASE + "/" + state.datasetId + "/expr_" +
      String(chunk).padStart(3, "0") + ".bin", Uint8Array);
    state.exprChunks[key] = p;
    return p.then(function (u8) {
      var off = (geneIdx - chunk * meta.chunkSize) * meta.nSpots;
      var vmax = meta.exprMax[geneIdx];
      var out = new Float32Array(meta.nSpots);
      for (var i = 0; i < meta.nSpots; i++) out[i] = (u8[off + i] / 255) * vmax;
      return out;
    });
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
    var want = state.token;
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

  function renderPlots() {
    var meta = state.meta;
    var gi = state.geneIdx;
    var gene = meta.genes[gi];
    var whole = state.cellType === "whole";

    plots.className = "dbx-plots" + (whole ? "" : " dbx-two");
    var propPanel = el("#dbx-prop-title").closest(".dbx-panel");
    propPanel.style.display = whole ? "none" : "";

    el("#dbx-expr-title").textContent = gene + " expression";
    el("#dbx-expr-sub").textContent =
      "Normalised count per spot · " + meta.label;

    var want = ++state.token;
    geneValues(gi).then(function (vals) {
      if (want !== state.token) return;
      draw("dbx-expr", vals, meta.exprMax[gi] || 1, "blue", "Normalised count");

      if (!whole) {
        var ctIdx = meta.cellTypes.indexOf(state.cellType);
        var nCT = meta.cellTypes.length - 1;          // cellProp has no "whole"
        var col = ctIdx - 1;
        var props = new Float32Array(meta.nSpots);
        var pmax = 0;
        for (var i = 0; i < meta.nSpots; i++) {
          props[i] = state.cellProp[i * nCT + col];
          if (props[i] > pmax) pmax = props[i];
        }
        el("#dbx-prop-title").textContent = state.cellType + " proportion";
        el("#dbx-prop-sub").textContent =
          "Estimated proportion per spot · " + meta.label;
        draw("dbx-prop", props, pmax || 1, "orange", "Proportion");
      }
      setStatus("");
    }).catch(function (e) {
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
        return '<div role="option" data-i="' + i + '" aria-selected="false">' + g + "</div>";
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

  /* Prefix matches first - typing "GFA" should surface GFAP before AGFAP-like
     substring hits - then substring matches, capped so the list stays usable. */
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
  new MutationObserver(function () {
    if (state.meta && state.geneIdx >= 0) renderPlots();
  }).observe(document.documentElement, { attributes: true, attributeFilter: ["data-bs-theme"] });

  var resizeTimer;
  window.addEventListener("resize", function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () {
      if (state.meta && state.geneIdx >= 0) renderPlots();
    }, 150);
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
        return '<option value="' + c + '">' + c + "</option>";
      }).join("");
      if (state.meta.cellTypes.indexOf(state.cellType) < 0) state.cellType = "whole";
      ctSel.value = state.cellType;

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

  Promise.all([
    fetchJson(BASE + "/manifest.json"),
    fetchJson(BASE + "/pops_genes.json")
  ]).then(function (res) {
    state.manifest = res[0];
    state.popsGenes = res[1];

    el("#dbx-dataset").innerHTML = state.manifest.datasets.map(function (d) {
      return '<option value="' + d.id + '">' + d.label + "</option>";
    }).join("");

    el("#dbx-disease").innerHTML = '<option value="-1">None</option>' +
      state.manifest.diseases.map(function (d, i) {
        return '<option value="' + i + '">' + d + "</option>";
      }).join("");

    return loadDataset(state.manifest.datasets[0].id);
  }).catch(function (e) {
    setStatus("Could not load the data browser: " + e.message, true);
  });
})();
