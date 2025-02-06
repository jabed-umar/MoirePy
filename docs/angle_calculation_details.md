# Angle Calculation Details

This will have information about given some certain a and b how the angle is calculated, why are the calculations such? how to calculate the moire lattice vector given a and b... why? illustrate with diagram and examples....


given an a and a b how the angle is calculated @jabed write. Here is a javascript code for reference (note that this loops through all start and end... you just need to mathematically write how given an a and a b, the theta can be calculated):

**Inputs:**  
- \( a, b \) (integer coefficients)  
- **Layer 1 Lattice Vectors:** \( \mathbf{A_1}, \mathbf{B_1} \)  
- **Layer 2 Lattice Vectors:** \( \mathbf{A_2}, \mathbf{B_2} \)  

**Vectors:**  
\[
\mathbf{V_1} = a \mathbf{A_1} + b \mathbf{B_1}
\]
\[
\mathbf{V_2} = b \mathbf{A_2} + a \mathbf{B_2}
\]

**Angle \( \theta \) Between \( \mathbf{V_1} \) and \( \mathbf{V_2} \):**  
\[
\theta = \cos^{-1} \left( \frac{\mathbf{V_1} \cdot \mathbf{V_2}}{||\mathbf{V_1}|| \cdot ||\mathbf{V_2}||} \right)
\]  
(Ensure \( \gcd(a, b) = 1 \) for valid results.)


```javascript
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
```