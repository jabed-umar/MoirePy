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
	// print the Number of results
	console.log("Number of results: ", results.length);
	displayResults(results);
}

	
function find_values(start, end, layer1Vectors, layer2Vectors) {
    const [a1x, a1y, b1x, b1y] = layer1Vectors;
    const [a2x, a2y, b2x, b2y] = layer2Vectors;
    const dot = (v1, v2) => v1[0] * v2[0] + v1[1] * v2[1];
    const norm = (v) => Math.sqrt(v[0] * v[0] + v[1] * v[1]);
    const results = [];
    for (let a = start; a <= end; a++) {
        for (let b = start; b <= end; b++) {
			// if a < b or b < 1 then skip
			if (a >= b || a < 1) continue;
            const one = [a * a1x + b * b1x, a * a1y + b * b1y];
            const two = [b * a2x + a * b2x, b * a2y + a * b2y];
            const c = dot(one, two) / (norm(one) * norm(two));
            const thetaRad = Math.acos(c);
            const thetaDeg = (thetaRad * 180) / Math.PI;
			if (gcd(a, b) !== 1) continue;
            results.push([thetaDeg.toFixed(8), thetaRad.toFixed(8), a, b]);
        }
    }
    results.sort((x, y) => x[0] - y[0]);
    return results;
}


function displayResults(results) {
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