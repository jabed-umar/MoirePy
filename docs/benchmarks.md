# MoirePy Benchmarks

This page presents interactive benchmark results for MoirePy across multiple machines, using the same benchmark protocol and plotting format so comparisons stay fair and readable. You can inspect **full-pipeline (lattice points + connections + hamiltonian)** time and **Hamiltonian-only** time as lattice size grows. The table below each plot summarizes platform metadata and headline timings at the largest sampled size. As more systems are added, this page becomes a living performance map of real-world hardware. If you want to add your machine, see the [Benchmark Contribution Guide](benchmark_contribution.md).

??? info "How to use the interactive plots"
    - Drag to zoom
    - Double-click to reset zoom
    - Toggle systems using checkboxes
    - Hover to inspect values

<div id="bench-status" style="margin: 0.5rem 0 1rem 0;"></div>

<div id="bench-controls" style="display:flex; flex-wrap:wrap; gap:0.75rem; margin:0.75rem 0 1.25rem 0;"></div>

## Generate Lattice Points, Connections, and Hamiltonian

<div id="chart-exp1" style="width:100%; min-height:360px;"></div>

## Subsequent Hamiltonian Generation after connections are made

<div id="chart-exp2" style="width:100%; min-height:360px;"></div>

## Systems

<table id="bench-systems">
  <thead>
    <tr>
      <th>Label</th>
      <th>CPU</th>
      <th>Platform</th>
      <th>N</th>
      <th>Repeat</th>
      <th>Max points</th>
      <th>Exp1 @ max (ms)</th>
      <th>Exp2 @ max (ms)</th>
    </tr>
  </thead>
  <tbody></tbody>
</table>

<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/uplot@1.6.31/dist/uPlot.min.css">
<script src="https://cdn.jsdelivr.net/npm/uplot@1.6.31/dist/uPlot.iife.min.js"></script>

<script>
const DATASETS = [
  { label: "i5 Gen12", path: "../benchmark/i5_gen12/benchmark_results.json", color: "#1472e6", enabled: true },
  { label: "i3 Gen4", path: "../benchmark/i3_gen4/benchmark_results.json", color: "#e95d2a", enabled: false },
  { label: "Raspberry Pi 5", path: "../benchmark/pi5/benchmark_results.json", color: "#0f9d58", enabled: false },
];

const statusEl = document.getElementById("bench-status");
const controlsEl = document.getElementById("bench-controls");
const tableBody = document.querySelector("#bench-systems tbody");
const exp1Host = document.getElementById("chart-exp1");
const exp2Host = document.getElementById("chart-exp2");

let exp1Plot = null;
let exp2Plot = null;
let loaded = [];

const fmt = (v, digits = 3) => (v == null ? "-" : Number(v).toFixed(digits));

function setStatus(text, isError = false) {
  statusEl.textContent = text;
  statusEl.style.color = isError ? "#d83a3a" : "inherit";
}

async function loadAll() {
  setStatus("Loading benchmark datasets...");

  const results = await Promise.all(
    DATASETS.map(async (ds) => {
      try {
        const res = await fetch(ds.path, { cache: "no-store" });
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const json = await res.json();
        return { ...ds, json, ok: true };
      } catch (err) {
        return { ...ds, ok: false, error: String(err) };
      }
    })
  );

  const failed = results.filter((r) => !r.ok);
  loaded = results.filter((r) => r.ok);

  if (loaded.length === 0) {
    setStatus("Failed to load all benchmark JSON files.", true);
    return;
  }

  if (failed.length > 0) {
    setStatus(`Loaded ${loaded.length} dataset(s). Failed: ${failed.map((f) => f.label).join(", ")}`, true);
  } else {
    setStatus(`Loaded ${loaded.length} benchmark dataset(s).`);
  }

  renderControls();
  renderTable();
  renderCharts();
}

function renderControls() {
  controlsEl.innerHTML = "";

  loaded.forEach((ds, i) => {
    const id = `bench-toggle-${i}`;
    const wrapper = document.createElement("label");
    wrapper.setAttribute("for", id);
    wrapper.style.display = "inline-flex";
    wrapper.style.alignItems = "center";
    wrapper.style.gap = "0.4rem";
    wrapper.style.padding = "0.35rem 0.55rem";
    wrapper.style.border = "1px solid #d7d7d7";
    wrapper.style.borderRadius = "0.45rem";
    wrapper.style.cursor = "pointer";

    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.id = id;
    cb.checked = Boolean(ds.enabled);
    cb.dataset.index = String(i);
    cb.addEventListener("change", renderCharts);

    const swatch = document.createElement("span");
    swatch.style.display = "inline-block";
    swatch.style.width = "0.7rem";
    swatch.style.height = "0.7rem";
    swatch.style.borderRadius = "999px";
    swatch.style.background = ds.color;

    const text = document.createElement("span");
    text.textContent = ds.label;

    wrapper.appendChild(cb);
    wrapper.appendChild(swatch);
    wrapper.appendChild(text);
    controlsEl.appendChild(wrapper);
  });
}

function getEnabled() {
  const boxes = controlsEl.querySelectorAll("input[type='checkbox']");
  return Array.from(boxes)
    .filter((b) => b.checked)
    .map((b) => loaded[Number(b.dataset.index)]);
}

function buildAlignedSeries(enabled, key) {
  const xSet = new Set();
  enabled.forEach((ds) => ds.json.results.forEach((r) => xSet.add(r.cells_times_2)));
  const xVals = Array.from(xSet).sort((a, b) => a - b);

  const ys = enabled.map((ds) => {
    const map = new Map(ds.json.results.map((r) => [r.cells_times_2, r[key]]));
    return xVals.map((x) => (map.has(x) ? map.get(x) : null));
  });

  return [xVals, ...ys];
}

function baseOpts(title, enabled, hostEl) {
  const width = Math.max(640, hostEl.clientWidth || 640);

  return {
    title,
    width,
    height: 360,
    scales: {
      x: { time: false },
      y: { auto: true },
    },
    axes: [
      { label: "Number of lattice points" },
      { label: "Time taken (ms)" },
    ],
    cursor: { drag: { x: true, y: true } },
    legend: { show: true },
    series: [
      { label: "points" },
      ...enabled.map((ds) => ({ label: ds.label, stroke: ds.color, width: 2 })),
    ],
  };
}

function renderCharts() {
  const enabled = getEnabled();

  if (enabled.length === 0) {
    if (exp1Plot) {
      exp1Plot.destroy();
      exp1Plot = null;
    }
    if (exp2Plot) {
      exp2Plot.destroy();
      exp2Plot = null;
    }
    exp1Host.innerHTML = "<p>Select at least one dataset.</p>";
    exp2Host.innerHTML = "<p>Select at least one dataset.</p>";
    return;
  }

  exp1Host.innerHTML = "";
  exp2Host.innerHTML = "";

  const dataExp1 = buildAlignedSeries(enabled, "experiment_1_ms");
  const dataExp2 = buildAlignedSeries(enabled, "experiment_2_ms");

  if (exp1Plot) exp1Plot.destroy();
  if (exp2Plot) exp2Plot.destroy();

  exp1Plot = new uPlot(baseOpts("Experiment 1: Full Pipeline", enabled, exp1Host), dataExp1, exp1Host);
  exp2Plot = new uPlot(baseOpts("Experiment 2: Hamiltonian Only", enabled, exp2Host), dataExp2, exp2Host);
}

function renderTable() {
  tableBody.innerHTML = "";

  loaded.forEach((ds) => {
    const config = ds.json.config || {};
    const results = ds.json.results || [];
    const last = results[results.length - 1] || {};

    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${ds.label}</td>
      <td>${config.cpu_name || "-"}</td>
      <td>${config.platform_system || "-"} ${config.platform_machine || ""}</td>
      <td>${config.n ?? "-"}</td>
      <td>${config.repeat ?? "-"}</td>
      <td>${last.cells_times_2 ?? "-"}</td>
      <td>${fmt(last.experiment_1_ms)}</td>
      <td>${fmt(last.experiment_2_ms)}</td>
    `;
    tableBody.appendChild(tr);
  });
}

window.addEventListener("resize", () => {
  if (!exp1Plot && !exp2Plot) return;
  renderCharts();
});

loadAll();
</script>
