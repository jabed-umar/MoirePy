<style>

.input-form {
    display: flex;
    flex-direction: column;
    gap: 10px;
    margin-bottom: 20px;
}

input, select, button {
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

table {
    width: 100%;
    border-collapse: collapse;
    margin-top: 20px;
}

th, td {
    padding: 8px;
    text-align: center;
    border: 1px solid #ddd;
}

thead {
    background-color: #f2f2f2;
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

</style>

# Angle Value Calculator

The `Moire` class takes two integers, `a` and `b`, to compute `θ` (the rotation angle of the upper lattice while keeping the lower lattice fixed) and other quantities like moiré lattice vectors. This tool helps you find suitable values of `a` and `b` for your desired angles.
*(Finding `a` and `b` for a given `θ` is inefficient, as not all rotations are possible for finite lattice sizes.)*

## Guidelines

- Currently, **both layers must be the same**. We do not support different layer types yet, but we plan to in the future.
- Results depend **only on lattice vectors**, not layer types. Identical lattice vectors (e.g., `HexagonalLayer` and `TriangularLayer`) will yield the same results.

## Choosing Start and End Values

- `Start` and `End` define the range for `a` and `b`. The tool computes `θ` for all valid pairs within this range.
- A **higher LCM** of `a` and `b` results in a **smaller θ**. **Larger values of `a` and `b`** produce even smaller angles, but this significantly increases the number of results. Expect longer processing times.

For details on how `θ` is calculated, see [Angle Calculation Details](angle_calculation_details.md).


## Calculator

<div class="input-form">
<div style="display: flex; gap: 10px; align-items: center;">
    <label for="start">Start:</label>
    <input type="number" id="start" value="0" style="width: 30vw;">

    <label for="end">End:</label>
    <input type="number" id="end" value="20" style="width: 30vw;">
</div>



 <!-- Layer 1 Selection -->
 <label for="layer1">Layer 1:</label>
 <select id="layer1" onchange="updateLatticeVectors('layer1')">
  <option value="Custom">Custom</option>
  <option value="SquareLayer">SquareLayer</option>
  <option value="RhombusLayer">RhombusLayer</option>
  <option value="TriangularLayer" selected>TriangularLayer</option>
  <option value="HexagonalLayer">HexagonalLayer</option>
  <option value="KagomeLayer">KagomeLayer</option>
 </select>

 <div id="layer1-vectors" class="vector-inputs hidden">
  <label>Lattice Vectors for Layer 1:</label>
  <div>lv1: <input type="text" inputmode="decimal" id="layer1-lv1x"> x + <input type="text" inputmode="decimal"
    id="layer1-lv1y"> y</div>
  <div>lv2: <input type="text" inputmode="decimal" id="layer1-lv2x"> x + <input type="text" inputmode="decimal"
    id="layer1-lv2y"> y</div>
 </div>

 <!-- Layer 2 Selection -->
 <label for="layer2">Layer 2:</label>
 <select id="layer2" onchange="updateLatticeVectors('layer2')">
  <option value="Custom">Custom</option>
  <option value="SquareLayer">SquareLayer</option>
  <option value="RhombusLayer">RhombusLayer</option>
  <option value="TriangularLayer" selected>TriangularLayer</option>
  <option value="HexagonalLayer">HexagonalLayer</option>
  <option value="KagomeLayer">KagomeLayer</option>
 </select>

 <div id="layer2-vectors" class="vector-inputs hidden">
  <label>Lattice Vectors for Layer 2:</label>
  <div>lv1: <input type="text" inputmode="decimal" id="layer2-lv1x"> x + <input type="text" inputmode="decimal"
    id="layer2-lv1y"> y</div>
  <div>lv2: <input type="text" inputmode="decimal" id="layer2-lv2x"> x + <input type="text" inputmode="decimal"
    id="layer2-lv2y"> y</div>
 </div>

 <button onclick="calculate()">CALCULATE</button>
</div>

<table id="results-table">
 <thead>
  <tr>
   <th>angle (deg)</th>
   <th>angle (rad)</th>
   <th>a</th>
   <th>b</th>
  </tr>
 </thead>
 <tbody id="results-body">
  <!-- Generated results will be displayed here -->
 </tbody>
</table>
</div>

<!-- <script src="assets/script_find_theta.js"></script> -->


<script>
let root3 = Math.sqrt(3);

const latticeDefaults = {
    HexagonalLayer:  [1, 0, 0.5, root3 / 2],
    SquareLayer:     [1, 0,   0,         1],
    RhombusLayer:    [1, 0, 0.5, root3 / 2],
    TriangularLayer: [1, 0, 0.5, root3 / 2],
    KagomeLayer:     [1, 0, 0.5, root3 / 2],
};

function updateLatticeVectors(layerId) {
    // only for custom vectors
    const layer = document.getElementById(layerId).value;
    const vectorContainer = document.getElementById(`${layerId}-vectors`);
    if (layer === "Custom") {
        vectorContainer.classList.remove("hidden");
    } else {
        vectorContainer.classList.add("hidden");
        const vectors = latticeDefaults[layer] || [1, 1, 1, 1]; // Placeholder if values are unknown
        document.getElementById(`${layerId}-lv1x`).value = vectors[0];
        document.getElementById(`${layerId}-lv1y`).value = vectors[1];
        document.getElementById(`${layerId}-lv2x`).value = vectors[2];
        document.getElementById(`${layerId}-lv2y`).value = vectors[3];
    }
}

function gcd(x, y) {
    if (y === 0) return x;
    else return gcd(y, x % y);
}

function calculate() {
    const start = parseInt(document.getElementById("start").value);
    const end = parseInt(document.getElementById("end").value);
    const layer1Vectors = [
        parseFloat(document.getElementById("layer1-lv1x").value),
        parseFloat(document.getElementById("layer1-lv1y").value),
        parseFloat(document.getElementById("layer1-lv2x").value),
        parseFloat(document.getElementById("layer1-lv2y").value)
    ];
    const layer2Vectors = [
        parseFloat(document.getElementById("layer2-lv1x").value),
        parseFloat(document.getElementById("layer2-lv1y").value),
        parseFloat(document.getElementById("layer2-lv2x").value),
        parseFloat(document.getElementById("layer2-lv2y").value)
    ];

    const results = find_values(start, end, layer1Vectors, layer2Vectors);
    console.log("Number of results: ", results.length);
    console.log(displayResults);
    displayResults_(results);
}

function find_values(start, end, layer1Vectors, layer2Vectors) {
    const [a1x, a1y, b1x, b1y] = layer1Vectors;
    const [a2x, a2y, b2x, b2y] = layer2Vectors;
    const dot = (v1, v2) => v1[0] * v2[0] + v1[1] * v2[1];
    const norm = (v) => Math.sqrt(v[0] * v[0] + v[1] * v[1]);
    const results = [];
    for (let a = start; a <= end; a++) {
        for (let b = start; b <= end; b++) {
            if (a >= b || a < 1) continue;  // checks
            const one = [a * a1x + b * b1x, a * a1y + b * b1y];
            const two = [b * a2x + a * b2x, b * a2y + a * b2y];
            const c = dot(one, two) / (norm(one) * norm(two));
            const thetaRad = Math.acos(c);
            const thetaDeg = (thetaRad * 180) / Math.PI;
            if (gcd(a, b) !== 1) continue;  // checks
            results.push([thetaDeg.toFixed(8), thetaRad.toFixed(8), a, b]);
        }
    }
    results.sort((x, y) => x[0] - y[0]);
    return results;
}

function displayResults_(results) {
    console.log(results);
    const resultsBody = document.getElementById("results-body");
    resultsBody.innerHTML = ""; // Clear previous results
    results.forEach(tuple => {
        const row = document.createElement("tr");
        tuple.forEach(value => {
            const cell = document.createElement("td");
            cell.textContent = value;
            row.appendChild(cell);
        });
        resultsBody.appendChild(row);
    });
}

updateLatticeVectors('layer1');
updateLatticeVectors('layer2');

</script>