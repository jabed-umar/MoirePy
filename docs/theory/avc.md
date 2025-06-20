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

    table {
        width: 100%;
        border-collapse: collapse;
        margin-top: 20px;
    }

    th,
    td {
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
    .section{
        text-align: justify;
    }
</style>





# Angle Value Calculator

The **Moiré Angle Calculator** finds all possible commensurate angles between two stacked lattices by trimming both layers within a given radius and checking for periodic overlaps. For each valid angle, it returns not just the angle itself, but also the corresponding lattice coordinates: `ll1, ll2` (from the lower lattice) and `ul1, ul2` (from the upper one). Since most angles are irrational and can't be precisely represented, the rotation is defined using these coordinate pairs instead — they mark the overlapping points between the two lattices. Once you've selected an angle from the list, just copy the `ll1, ll2, ul1, ul2` values into your code — the system will figure out the rotation angle from that. For the logic behind how angles are identified, see the [Angle Calculation Process](angle_calculation_process.md).


## Guidelines

<details>
    <summary>Click to Expand</summary>

    <ul>
        <li>
            A <strong>radius</strong> must be specified to define the extent of the circular region centered at the origin. This radius will be used to truncate both the upper and lower lattices.
        </li>

        <li>
            <strong>Currently Supported Systems</strong>:
            <ul>
                <li>Triangular on triangular lattices</li>
                <li>Square on square lattices</li>
                <li>
                    A <strong>custom mode</strong> is available, allowing the input of arbitrary lattice vectors for each layer. Please note that this mode is experimental and its <strong>reliability is NOT guaranteed</strong>. Potential issues may include:
                    <ul>
                        <li>Erroneous or nonsensical output.</li>
                        <li>Unresponsiveness or infinite loops.</li>
                        <li>Unexpected program behavior.</li>
                    </ul>
                    However, if you believe that your specific configuration (as detailed in <a href="../angle_calculation_process">Angle Calculation Process</a>) <em>should</em> yield meaningful results, then following these patterns <em>might</em> be considered:
                    <ol>
                        <li>Both lattice angles are exact divisors of 360°.</li>
                        <li>One of them should exactly divide the other.</li>
                    </ol>
                    Even under these conditions, there is <strong>no assurance</strong> that the calculated angles will be accurate or relevant. Use this feature at your own discretion.
                </li>
            </ul>
        </li>

        <li>A <strong>larger radius</strong> will encompass more lattice points, potentially leading to more precise calculations and the detection of smaller angles, but will increase computation time.</li>
        <li>A <strong>smaller radius</strong> will yield faster results, but may only reveal larger angle values.</li>

    </ul>
</details>




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

    <div style="display: flex; gap: 10px; align-items: center;">
        <label for="precision">Precision:</label>
        <input type="range" id="precision" min="2" max="12" value="6" step="1" style="width: 30vw;">
        <span id="precision-value">6</span>
    </div>   

    <button onclick="calculate()">CALCULATE</button>
</div>



**Note:** The last column (number of points per unit cell in the moire lattice) has been calculated assuming only one point per unit cell. If you are using lattices which have multiple points per unit cell like hexagonal (2) or kagome (3), multiply this value by the number of points per unit cell in your lattice to get the actual number of points in the moire lattice.

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
            <th>points/cell</th>
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
        HexagonalLayer: [1, 0, 0.5, root3 / 2],
        SquareLayer: [1, 0, 0, 1],
        RhombusLayer: [1, 0, 0.5, root3 / 2],
        TriangularLayer: [1, 0, 0.5, root3 / 2],
        KagomeLayer: [1, 0, 0.5, root3 / 2],
    };

    document.getElementById("precision").addEventListener("input", function() {
        document.getElementById("precision-value").textContent = this.value;
    });

    function updateLatticeVectors() {
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

    function calculate() {
        const radius = parseInt(document.getElementById("radius").value);

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

        const precision = parseInt(document.getElementById("precision").value);

        // console.log(radius, layer1Vectors, layer2Vectors);

        const results = find_values(radius, layer1Vectors, layer2Vectors, tol=precision);

        console.log(results);
        console.log("Number of results:", results.length);
        displayResults_(results);
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

    function find_values(radius, layer1Vectors, layer2Vectors, tol = 6) {

        const [a1x, a1y, b1x, b1y] = layer1Vectors;
        const [a2x, a2y, b2x, b2y] = layer2Vectors;

        if (JSON.stringify(layer1Vectors) !== JSON.stringify(layer2Vectors)) {
            alert("Warning: Vectors are not identical! Results might be inaccurate or meaningless.");
        }

        const lv1 = [a1x, a1y];
        const lv2 = [b1x, b1y];

        const lattice1 = generate_lattice_points(lv1, lv2, radius);
        const lattice2 = generate_lattice_points(lv1, lv2, radius);

        const [dict1, dist_set1] = process_lattice(lattice1, tol);
        const [dict2, dist_set2] = process_lattice(lattice2, tol);

        const common_dists = [...dist_set1].filter(d => dist_set2.has(d)).sort((a, b) => a - b).slice(1);

        const angle_dict = {};
        const lattice_angle = angle_from_x(lv2) - angle_from_x(lv1);

        const isValidTheta = (theta) => theta > 0 && theta < lattice_angle;

        for (const d of common_dists) {
            // console.log(d)
            const pts1 = Object.values(dict1.get(d)).filter(p => isValidTheta(angle_from_x(p)));
            const pts2 = Object.values(dict2.get(d)).filter(p => isValidTheta(angle_from_x(p)));

            for (const p1 of pts1) {
                const theta1 = parseFloat(angle_from_x(p1).toFixed(tol));

                for (const p2 of pts2) {
                    const theta2 = parseFloat(angle_from_x(p2));
                    const angle = parseFloat((theta2 - theta1));
                    // use cos theta square between p1 and p2 as uid
                    const uid = angleId(p1, p2);

                    if (
                        theta2 <= theta1 ||
                        angle < Math.pow(10, -tol) ||
                        uid in angle_dict
                    ) continue;

                    angle_dict[uid] = [p1, p2, angle];
                }
            }
        }

        const results = Object.entries(angle_dict)
        .sort(([, a], [, b]) => parseFloat(a[2]) - parseFloat(b[2]))  // ascending by angle
        .map(([k, [p1, p2, angle]]) => {
            const thetaRad = (parseFloat(angle) * Math.PI) / 180;
            const thetaDeg = parseFloat(angle);
            const [i1, j1] = calc_indices(p1, lv1, lv2);
            const [i2, j2] = calc_indices(p2, lv1, lv2);
            const num_pts =   2*(p1[0] * p1[0] + p1[1] * p1[1]) * 1;  // 1 for one point per unit cell
            return [thetaDeg.toFixed(tol), thetaRad.toFixed(tol), i2, j2, i1, j1, num_pts];
        });

        return results;
    }

    function displayResults_(results) {
        const resultsBody = document.getElementById("results-body");
        resultsBody.innerHTML = ""; // Clear previous results
        results.forEach((tuple, index) => {
            const row = document.createElement("tr");
            const cell = document.createElement("td");
            cell.textContent = index + 1;  // add the index
            row.appendChild(cell);
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