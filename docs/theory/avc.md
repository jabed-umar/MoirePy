<style>
    .input-form {
        display: flex;
        flex-direction: column;
        gap: 10px;
        margin-bottom: 20px;
    }

    input,
    select,
    button {
        font-size: 16px;
        padding: 5px;
    }

    button {
        cursor: pointer;
        background-color: rgb(41, 128, 185);
        color: white;
        border: none;
        padding: 10px;
    }

    /* ===== Custom table styles for avc.md only ===== */

    /* Override Material theme's min-width on th (it sets min-width: 5rem via .md-typeset table:not([class]) th) */
    .md-typeset #results-table th,
    .md-typeset #results-table td {
        min-width: 0;
        padding: 8px 6px;
        text-align: center;
        word-break: break-word;
    }

    .md-typeset #results-table {
        width: 30rem;
        border-collapse: collapse;
        table-layout: fixed;
        margin-top: 20px;
    }

    /* Column widths: # | deg | rad | ll1 | ll2 | ul1 | ul2 | points/cell */
    .md-typeset #results-table th:nth-child(1),
    .md-typeset #results-table td:nth-child(1) { width: 5%; }
    .md-typeset #results-table th:nth-child(2),
    .md-typeset #results-table td:nth-child(2) { width: 10%; }
    .md-typeset #results-table th:nth-child(3),
    .md-typeset #results-table td:nth-child(3) { width: 10%; }
    .md-typeset #results-table th:nth-child(4),
    .md-typeset #results-table td:nth-child(4) { width: 9%; }
    .md-typeset #results-table th:nth-child(5),
    .md-typeset #results-table td:nth-child(5) { width: 9%; }
    .md-typeset #results-table th:nth-child(6),
    .md-typeset #results-table td:nth-child(6) { width: 9%; }
    .md-typeset #results-table th:nth-child(7),
    .md-typeset #results-table td:nth-child(7) { width: 9%; }
    .md-typeset #results-table th:nth-child(8),
    .md-typeset #results-table td:nth-child(8) { width: 10%; }

    /* Light mode colors */
    [data-md-color-scheme="default"] .md-typeset #results-table th {
        background: #f5f5f5;
        color: #222;
    }
    [data-md-color-scheme="default"] .md-typeset #results-table td {
        border-color: #ddd;
    }
    [data-md-color-scheme="default"] .md-typeset #results-table tr:nth-child(even) {
        background: #fafafa;
    }
    [data-md-color-scheme="default"] .md-typeset #results-table tr:hover {
        background: #eaeaea;
    }

    /* Dark mode colors */
    [data-md-color-scheme="slate"] .md-typeset #results-table th {
        background: #222;
        color: #f5f5f5;
    }
    [data-md-color-scheme="slate"] .md-typeset #results-table td {
        border-color: #444;
    }
    [data-md-color-scheme="slate"] .md-typeset #results-table tr:nth-child(even) {
        background: #2a2a2a;
    }
    [data-md-color-scheme="slate"] .md-typeset #results-table tr:hover {
        background: #333;
    }

    .hidden {
        display: none;
    }

    .vector-inputs input {
        width: 80px;
        margin: 5px;
    }

    #end, #start {
        border: 1px solid #000;
        border-radius: 5px;
        padding: 8px;
        width: 100px;
        font-size: 16px;
    }
    .section{
        text-align: justify;
    }

    #avc-context-menu {
        position: fixed;
        z-index: 10000;
        min-width: 170px;
        border-radius: 8px;
        overflow: hidden;
        box-shadow: 0 10px 24px rgba(0, 0, 0, 0.25);
    }

    #avc-context-menu button {
        width: 100%;
        text-align: left;
        border-radius: 0;
        background: transparent;
        color: inherit;
        padding: 10px 12px;
    }

    #avc-context-menu button:hover {
        background: rgba(127, 127, 127, 0.2);
    }

    .md-typeset #results-body tr {
        cursor: context-menu;
    }

    [data-md-color-scheme="default"] #avc-context-menu {
        background: #ffffff;
        color: #222;
        border: 1px solid #ddd;
    }

    [data-md-color-scheme="slate"] #avc-context-menu {
        background: #1f1f1f;
        color: #f5f5f5;
        border: 1px solid #444;
    }
</style>





# Angle Value Calculator

Most physically relevant twist configurations are irrational, so typing a rounded angle directly into code is usually unstable and hard to reproduce.

This calculator instead returns exact integer tuples `ll1`, `ll2`, `ul1`, `ul2` that define the commensurate match between the lower and upper layers. Those integers can be used directly in MoirePy and provide deterministic geometry reconstruction.

Internally, both layers are truncated inside a chosen radius and overlap candidates are filtered using the [Angle Calculation Process](angle_calculation_process.md).  
Each result row gives:

- angle in degrees and radians
- integer tuple (`ll1`, `ll2`, `ul1`, `ul2`)
- estimated supercell size (`cells`)


??? warning "Some Guidelines"
    A **radius** must be specified to define the extent of the circular region centred at the origin. This value is used to truncate *both* the upper and lower lattices. **Larger radius** Includes more lattice points, potentially giving more precise calculations and revealing smaller angles, **but** increases computation time. **Smaller radius** Produces faster results, yet may detect only larger angle values.

    **Currently supported systems:**

    - Both layers **60$^\circ$** (Triangular, Hexagonal and Kagome lattices)
    - Both layers **90$^\circ$** (Square lattice)
    - **Custom mode** (*experimental*)
        
        *Possible problems*  
        - Erroneous or nonsensical output  
        - Unresponsiveness or infinite loops  
        - Unexpected program behaviour  

        Proceed only if you believe the results yielded by the [Angle Calculation Process](angle_calculation_process.md) are correct and meaningful for your custom lattice vectors.

    **Larger Value of points per cell** More points per cell in the lattice means larger hamiltonian matrix, which can lead to longer computation times and higher memory usage. Especially if planning to invert the hamiltonian matrix and find eigenvalues and eigenvectors later.

    **Number of Points** has been calculated assuming only one point per unit cell. If you are using lattices which have multiple points per unit cell like hexagonal (2) or kagome (3), multiply this value by the number of points per unit cell in your lattice to get the actual number of points in the Moiré lattice.

    **Precision** has **NO effect on the calculation**. It only affects the output format of the results. (In the calculation, we stick to integer values to avoid floating point errors).




## Calculator

<div class="input-form">
    <div style="display: flex; gap: 10px; align-items: center;">
        <label for="radius">Radius:</label>
        <input type="number" id="radius" value="20" style="width: 30vw;" step="any" min="0">
    </div>

    <label for="latticeType">Lattice Type:</label>
    <select id="latticeType" onchange="updateLatticeVectors()">
        <option value="TriangularLayer">both layers 60 degree (Triangular, Hexagonal and Kagome lattice)</option>
        <option value="SquareLayer">both layers 90 degree (Square lattice)</option>
        <option value="Custom">Custom (experimental)</option>
    </select>

    <div id="custom-vectors" class="vector-inputs hidden">
        <label>Lattice Vectors for Layer 1:</label>
        <div>
            lv1: <input type="text" inputmode="decimal" id="layer1-lv1x"> x +
            <input type="text" inputmode="decimal" id="layer1-lv1y"> y
        </div>
        <div>
            lv2: <input type="text" inputmode="decimal" id="layer1-lv2x"> x +
            <input type="text" inputmode="decimal" id="layer1-lv2y"> y
        </div>

        <label>Lattice Vectors for Layer 2:</label>
        <div>
            lv1: <input type="text" inputmode="decimal" id="layer2-lv1x"> x +
            <input type="text" inputmode="decimal" id="layer2-lv1y"> y
        </div>
        <div>
            lv2: <input type="text" inputmode="decimal" id="layer2-lv2x"> x +
            <input type="text" inputmode="decimal" id="layer2-lv2y"> y
        </div>
        <div style="color: red; border: 1px solid red; padding: 5px; text-align: center;"><b>WARNING:</b> RESULTS MIGHT BE MEANINGLESS OR INACCURATE</div>
    </div>

    <div style="display: flex; gap: 10px; align-items: center; flex-wrap: wrap;">
        <label>Angle range (degrees):</label>
        <div style="display: flex; gap: 6px; align-items: center;">
            <label for="angle-start">Start:</label>
            <input type="number" id="angle-start" value="0" min="0" max="360" step="any" style="width: 80px;">
        </div>
        <div style="display: flex; gap: 6px; align-items: center;">
            <label for="angle-end">End:</label>
            <input type="number" id="angle-end" value="60" min="0" max="360" step="any" style="width: 80px;">
        </div>
    </div>

    <div style="display: flex; gap: 10px; align-items: center;">
        <label for="precision">Precision (10^-x):</label>
        <input type="range" id="precision" min="2" max="12" value="6" step="1" style="width: 30vw;">
        <span id="precision-value">6</span>
    </div>

    <button onclick="calculate()">CALCULATE</button>
</div>


??? note

    The last column (**cells**) reports the total number of small cells in the moiré supercell including **both upper and lower layers**.

    If your lattice has multiple basis points per unit cell, scale accordingly:
    <ul>
        <li>Hexagonal (2 basis points): <code>actual_points = cells * 2</code></li>
        <li>Kagome (3 basis points): <code>actual_points = cells * 3</code></li>
    </ul>

You can right click (or long-touch) to copy rows as python code.

<table id="results-table">
    <thead>
        <tr>
            <th></th>
            <th>angle (deg)</th>
            <th>angle (rad)</th>
            <th>ll1</th>
            <th>ll2</th>
            <th>ul1</th>
            <th>ul2</th>
            <th title="Estimated supercell size">cells</th>
        </tr>
    </thead>
    <tbody id="results-body">
        <!-- Generated results will be displayed here -->
    </tbody>
</table>
</div>

<div id="avc-context-menu" class="hidden" role="menu" aria-label="Row actions">
    <button id="avc-copy-preset" type="button" role="menuitem" title="Copy ll1, ll2, ul1, ul2 tuple">Copy Row as python</button>
    <div style="height: 1px; background: rgba(127,127,127,0.2); margin: 4px 0;"></div>
    <button id="avc-download-csv" type="button" role="menuitem" title="Export entire table">Download table as csv</button>
</div>

<!-- <script src="assets/script_find_theta.js"></script> -->






<script type="module">

    // WASM engine — load once, non-blocking.
    // Path is relative to the rendered page (/theory/avc/ → /theory/assets/…)
    import init, { find_values as wasm_find_values }
        from '../avc-wasm/pkg/avc_wasm.js';

    let wasmReady = false;
    init().then(() => { wasmReady = true; });

    let root3 = Math.sqrt(3);

    const latticeDefaults = {
        HexagonalLayer: [1, 0, 0.5, root3 / 2],
        SquareLayer: [1, 0, 0, 1],
        RhombusLayer: [1, 0, 0.5, root3 / 2],
        TriangularLayer: [1, 0, 0.5, root3 / 2],
        KagomeLayer: [1, 0, 0.5, root3 / 2],
    };

    const contextMenu = document.getElementById("avc-context-menu");
    const copyPresetButton = document.getElementById("avc-copy-preset");
    const resultsBody = document.getElementById("results-body");
    let contextMenuRow = null;

    function closeContextMenu() {
        contextMenu.classList.add("hidden");
        contextMenuRow = null;
    }

    function getPresetStringFromRow(row) {
        const cells = row.querySelectorAll("td");
        if (cells.length < 8) return null;

        const ll1 = cells[3].textContent.trim();
        const ll2 = cells[4].textContent.trim();
        const ul1 = cells[5].textContent.trim();
        const ul2 = cells[6].textContent.trim();

        return `ll1=${ll1}, ll2=${ll2}, ul1=${ul1}, ul2=${ul2},`;
    }

    async function copyTextToClipboard(text) {
        try {
            await navigator.clipboard.writeText(text);
            return true;
        } catch (err) {
            const fallback = document.createElement("textarea");
            fallback.value = text;
            fallback.setAttribute("readonly", "");
            fallback.style.position = "absolute";
            fallback.style.left = "-9999px";
            document.body.appendChild(fallback);
            fallback.select();
            const copied = document.execCommand("copy");
            document.body.removeChild(fallback);
            return copied;
        }
    }

    resultsBody.addEventListener("contextmenu", (event) => {
        const row = event.target.closest("tr");
        if (!row) return;

        event.preventDefault();
        contextMenuRow = row;
        contextMenu.classList.remove("hidden");

        const menuWidth = 170;
        const menuHeight = 44;
        const padding = 8;
        const viewportWidth = window.innerWidth;
        const viewportHeight = window.innerHeight;

        let left = event.clientX;
        let top = event.clientY;

        if (left + menuWidth + padding > viewportWidth) {
            left = viewportWidth - menuWidth - padding;
        }
        if (top + menuHeight + padding > viewportHeight) {
            top = viewportHeight - menuHeight - padding;
        }

        contextMenu.style.left = `${Math.max(left, padding)}px`;
        contextMenu.style.top = `${Math.max(top, padding)}px`;
    });

    copyPresetButton.addEventListener("click", async () => {
        if (!contextMenuRow) return;
        const preset = getPresetStringFromRow(contextMenuRow);
        if (!preset) return;

        const success = await copyTextToClipboard(preset);
        closeContextMenu();

        if (success) {
            const original = copyPresetButton.textContent;
            copyPresetButton.textContent = "Copied!";
            setTimeout(() => {
                copyPresetButton.textContent = original;
            }, 900);
        } else {
            alert("Copy failed. Please copy the values manually.");
        }
    });

    const downloadCsvButton = document.getElementById("avc-download-csv");
    downloadCsvButton.addEventListener("click", () => {
        closeContextMenu();
        
        const table = document.getElementById("results-table");
        if (!table) return;

        let csvString = "";
        const rows = table.querySelectorAll("tr");
        
        for (const row of rows) {
            const cols = row.querySelectorAll("th, td");
            const rowData = [];
            for (const col of cols) {
                // Escape quotes and wrap in quotes if there's a comma
                let text = col.textContent.trim();
                if (text.includes(",") || text.includes('"')) {
                    text = `"${text.replace(/"/g, '""')}"`;
                }
                rowData.push(text);
            }
            csvString += rowData.join(",") + "\n";
        }

        const blob = new Blob([csvString], { type: "text/csv;charset=utf-8;" });
        const link = document.createElement("a");
        if (link.download !== undefined) {
            const url = URL.createObjectURL(blob);
            link.setAttribute("href", url);
            link.setAttribute("download", "commensurate_angles.csv");
            link.style.visibility = "hidden";
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            URL.revokeObjectURL(url);
        }
    });

    document.addEventListener("click", (event) => {
        if (!event.target.closest("#avc-context-menu")) {
            closeContextMenu();
        }
    });

    document.addEventListener("keydown", (event) => {
        if (event.key === "Escape") {
            closeContextMenu();
        }
    });

    window.addEventListener("scroll", closeContextMenu, true);
    window.addEventListener("resize", closeContextMenu);

    document.getElementById("precision").addEventListener("input", function() {
        document.getElementById("precision-value").textContent = this.value;
    });

    window.updateLatticeVectors = function() {
        const type = document.getElementById("latticeType").value;
        const customDiv = document.getElementById("custom-vectors");

        if (type === "Custom") {
            customDiv.classList.remove("hidden");
        } else {
            customDiv.classList.add("hidden");

            const vec = latticeDefaults[type] || [1, 1, 1, 1];

            // Set both layers with the same vectors
            document.getElementById("layer1-lv1x").value = vec[0];
            document.getElementById("layer1-lv1y").value = vec[1];
            document.getElementById("layer1-lv2x").value = vec[2];
            document.getElementById("layer1-lv2y").value = vec[3];

            document.getElementById("layer2-lv1x").value = vec[0];
            document.getElementById("layer2-lv1y").value = vec[1];
            document.getElementById("layer2-lv2x").value = vec[2];
            document.getElementById("layer2-lv2y").value = vec[3];
        }

        // Auto-set angle range based on lattice symmetry
        const angleEnd = (type === "SquareLayer") ? 90 : 60;
        document.getElementById("angle-start").value = 0;
        document.getElementById("angle-end").value = angleEnd;
    }

    function gcd(x, y) {
        if (y === 0) return x;
        else return gcd(y, x % y);
    }

    function angleId(p1, p2) {
        // Dot product
        const dot = p1[0] * p2[0] + p1[1] * p2[1];
        const dotSq = dot * dot;

        // Norms squared
        const norm1Sq = p1[0] ** 2 + p1[1] ** 2;
        const norm2Sq = p2[0] ** 2 + p2[1] ** 2;
        const denom = norm1Sq * norm2Sq;

        // Reduce the fraction dotSq / denom
        const commonDivisor = gcd(dotSq, denom);
        const num = dotSq / commonDivisor;
        const den = denom / commonDivisor;

        // Return as a string ID
        return `${num}/${den}`;
    }

    window.calculate = function() {
        if (!wasmReady) {
            alert("WASM engine is still loading — please try again in a moment.");
            return;
        }

        const radius = parseInt(document.getElementById("radius").value);

        const layer1Vectors = new Float64Array([
            parseFloat(document.getElementById("layer1-lv1x").value),
            parseFloat(document.getElementById("layer1-lv1y").value),
            parseFloat(document.getElementById("layer1-lv2x").value),
            parseFloat(document.getElementById("layer1-lv2y").value)
        ]);

        const layer2Vectors = new Float64Array([
            parseFloat(document.getElementById("layer2-lv1x").value),
            parseFloat(document.getElementById("layer2-lv1y").value),
            parseFloat(document.getElementById("layer2-lv2x").value),
            parseFloat(document.getElementById("layer2-lv2y").value)
        ]);

        const precision = parseInt(document.getElementById("precision").value);
        const angle_start = parseFloat(document.getElementById("angle-start").value);
        const angle_end = parseFloat(document.getElementById("angle-end").value);

        const t0 = performance.now();
        const flat = wasm_find_values(radius, layer1Vectors, layer2Vectors, precision, angle_start, angle_end);
        const t1 = performance.now();

        const N = flat.length / 7;
        console.log(`Found ${N} results in ${(t1 - t0).toFixed(2)} milliseconds.`);
        displayResults(flat, precision);
    }

    function calc_indices(p, lv1, lv2) {
        const [a, b] = lv1;
        const [c, d] = lv2;
        const [x, y] = p;
        const det = (a * d - b * c);
        const nx = (d * x - c * y) / det;
        const ny = (a * y - b * x) / det;

        if (Math.abs(Math.round(nx) - nx) > 1e-5 || Math.abs(Math.round(ny) - ny) > 1e-5) {
            throw new Error(`Calculation error for indices: ${nx}, ${ny}`);
        }

        return [Math.round(nx), Math.round(ny)];
    }

    function generate_lattice_points(lv1, lv2, radius) {
        const points = [];
        const maxGridSize = Math.floor(radius / Math.abs(lv2[1])) + 5;

        // console.log(radius, lv1, lv2, maxGridSize);

        for (let i = -maxGridSize; i <= maxGridSize; i++) {
            for (let j = -maxGridSize; j <= maxGridSize; j++) {
                // console.log(i, j);
                const point = [i * lv1[0] + j * lv2[0], i * lv1[1] + j * lv2[1]];
                const dist = Math.sqrt(point[0] ** 2 + point[1] ** 2);
                if (dist <= radius) points.push(point);
            }
        }

        return points;
    }

    function angle_from_x(p) {
        return Math.atan2(p[1], p[0]) * 180 / Math.PI;
    }

    function process_lattice(points, tol) {
        const distances = points.map(p => Math.hypot(p[0], p[1]));
        const distMap = new Map();

        // console.log(distMap);

        for (let i = 0; i < points.length; i++) {
            const d = parseFloat(distances[i].toFixed(tol));
            if (!distMap.has(d)) distMap.set(d, {});
            distMap.get(d)[i] = points[i];
        }

        return [distMap, new Set([...distMap.keys()])];
    }


    // flat: Float64Array with 7 values per row
    //   [ angle_deg, angle_rad, ll1, ll2, ul1, ul2, cells ]
    function displayResults(flat, precision) {
        resultsBody.innerHTML = "";
        const N = flat.length / 7;
        for (let i = 0; i < N; i++) {
            const row = document.createElement("tr");
            row.title = "Row actions available";
            const idxCell = document.createElement("td");
            idxCell.textContent = i + 1;
            row.appendChild(idxCell);
            for (let j = 0; j < 7; j++) {
                const cell = document.createElement("td");
                cell.textContent = j < 2
                    ? flat[i * 7 + j].toFixed(precision)
                    : Math.round(flat[i * 7 + j]);
                row.appendChild(cell);
            }
            resultsBody.appendChild(row);
        }
    }

    updateLatticeVectors();

</script>
