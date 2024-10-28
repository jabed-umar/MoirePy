let root3 = Math.sqrt(3);

function updateTolerance() {
	document.getElementById("tolerance-value").textContent = document.getElementById("tolerance").value;
}

const latticeDefaults = {
	HexagonalLayer: [1, 0, 0.5, root3 / 2],
	SquareLayer: [1, 0, 0, 1],
	RhombusLayer: [1, 0, 0.5, root3 / 2],
	TriangularLayer: [1, 0, 0.5, root3 / 2],
	KagomeLayer: [1, 0, 0.5, root3 / 2]
};

function updateLatticeVectors(layerId) {
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

function calculate() {
	const start = parseInt(document.getElementById("start").value);
	const end = parseInt(document.getElementById("end").value);
	const tolerance = parseFloat(document.getElementById("tolerance").value);

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

	const results = find_values(start, end, tolerance, layer1Vectors, layer2Vectors);
	// print the Number of results
	console.log("Number of results: ", results.length);
	displayResults(results);
}

function dotProduct(v1, v2) {
	return v1[0] * v2[0] + v1[1] * v2[1];
}

function norm(v) {
	return Math.sqrt(v[0] ** 2 + v[1] ** 2);
}

function angleBetweenVectors(v1, v2) {
	const cosTheta = dotProduct(v1, v2) / (norm(v1) * norm(v2));
	return Math.acos(cosTheta) * (180 / Math.PI);
}

function gcd_two_numbers(x, y) {
	if (y === 0) return x;
	else return gcd_two_numbers(y, x % y);
}


function check_abmn(a, b, m, n, a1, b1, a2, b2, tol) {
	// console.log("hello! ", Math.abs(Math.round(n) - n));
	// console.log(Math.abs(Math.round(n) - n), tol, n, Math.round(n));

	if (Math.abs(Math.round(n) - n) > 10 ** tol) return null;
	n = Math.round(n);

	// if (a === 0 || b === 0 || a === m || n <= 0) return null;
	if (a === 0 || b === 0 || a === m || n <= 0) return null;

	// if HCF of a and b is not 1, return null
	if (gcd_two_numbers(a, b) !== 1) return

	const vec1 = [a * a1[0] + b * b1[0], a * a1[1] + b * b1[1]];
	const vec2 = [m * a2[0] + n * b2[0], m * a2[1] + n * b2[1]];

	const thetadeg = angleBetweenVectors(vec1, vec2);
	const thetarad = thetadeg * Math.PI / 180;

	// if thetadeg is NaN or 0, return null
	if (isNaN(thetadeg) || thetadeg === 0) return null;

	return [thetadeg.toFixed(4), thetarad.toFixed(4), a, b, m, n];
}

function find_values(start, end, tolerance, layer1Vectors, layer2Vectors) {

	// console.log("start: ", start);
	// console.log("end: ", end);
	// console.log("tolerance: ", tolerance);
	// console.log("layer1Vectors: ", layer1Vectors);
	// console.log("layer2Vectors: ", layer2Vectors);

	const [a1x, a1y, b1x, b1y] = layer1Vectors;
	const [a2x, a2y, b2x, b2y] = layer2Vectors;
	const a1 = [a1x, a1y];
	const b1 = [b1x, b1y];
	const a2 = [a2x, a2y];
	const b2 = [b2x, b2y];

	const norm_a1 = norm(a1);
	const norm_b1 = norm(b1);
	const norm_a2 = norm(a2);
	const norm_b2 = norm(b2);

	const cos1 = dotProduct(a1, b1) / (norm_a1 * norm_b1);
	const cos2 = dotProduct(a2, b2) / (norm_a2 * norm_b2);
	const abs_n = norm_b2;
	// console.log("cos1: ", cos1);
	// console.log("cos2: ", cos2);
	// console.log("abs_n: ", abs_n);

	const results = [];

	for (let a = start + 1; a <= end; a++) {
		for (let b = start + 1; b <= a; b++) {
			for (let m = start + 1; m <= end; m++) {

				const aa = a * norm_a1;
				const bb = b * norm_b1;
				const mm = m * norm_a2;

				const sa_b = 2 * mm * cos2;
				const sa_c = mm ** 2 - (aa ** 2 + bb ** 2 + 2 * aa * bb * cos1);

				const discriminant = sa_b ** 2 - 4 * sa_c;
				if (discriminant < 0) continue;

				const sqrtDiscriminant = Math.sqrt(discriminant);

				let n = (-sa_b + sqrtDiscriminant) / (2 * abs_n);
				let result = check_abmn(a, b, m, n, a1, b1, a2, b2, tolerance);
				if (result) results.push(result);

				n = (-sa_b - sqrtDiscriminant) / (2 * abs_n);
				result = check_abmn(a, b, m, n, a1, b1, a2, b2, tolerance);
				if (result) results.push(result);
			}
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