<div style="text-align: justify;">
<h1>Angle Calculation Process</h1>

When two single layers of a 2D material are stacked on top of each other with a small misalignment, they produce a moiré pattern of much larger length scale than the periodicity of each layer. At a certain angle of twist, this results in a <em>commensurate moiré pattern</em>—structures where atoms from one layer exactly overlap atoms from the other. <br>

In this section, we describe how <strong>MoirePy</strong> determines the <em>commensurate rotation angles</em> between two lattices confined to a bounded region. The objective is to identify angles at which the two lattices align in a periodic fashion, resulting in well-ordered, physically observable moiré patterns.

</div>

<div style="text-align: justify;">
<div style="height: 10px;"></div>

<h2>Problem Statement</h2>

Let lattices \( A \) and \( B \) be two periodic point sets in two dimensions, each defined by their respective lattice vectors.

We are interested in the following geometric question:

Given a rotation angle \( \theta \), does there exist a point \( \mathbf{p} \in A \) and a point \( \mathbf{q} \in B \) such that
$$
\mathbf{p} = R(\theta)\mathbf{q}
$$

Here, \( R(\theta) \) denotes the standard two-dimensional rotation matrix:

$$
R(\theta) =
\begin{bmatrix}
\cos\theta & -\sin\theta \\
\sin\theta & \cos\theta
\end{bmatrix}
$$

To make the computation finite, we restrict our attention to a finite region by considering only lattice points within a circular region of radius \( r \). Let the truncated sets be:
\( A_r = A \cap \text{circle}(r) \) and \( B_r = B \cap \text{circle}(r) \).
<br>
Our objective is to determine the set of angles \( \theta \) for which there exists a pair of points \( \mathbf{p} \in A_r \) and \( \mathbf{q} \in B_r \) such that:

$$
\mathbf{p} = R(\theta)\mathbf{q}
$$

These angles correspond to potential commensurate alignments between the two lattices, giving rise to physically meaningful moiré patterns.
</div>

<div style="text-align: justify;">
<div style="height: 10px;"></div>

<h2>Traditional Diophantine Equation Approach</h2>

In commensurate moiré superlattices, there exist specific periodic points where atoms from the top and bottom layers align exactly. To analyze these alignments, we consider the primitive lattice vectors of the lower and upper layers as \( \vec{a}_1\), \(\vec{b}_1 \) and \( \vec{m}_1\), \(\vec{n}_1 \), respectively. The atomic positions in each layer are then given by:

\[
\vec{R}^1_{p,q} = a\vec{a}_1 + b\vec{b}_1 \quad \text{and} \quad
\vec{R}^2_{m,n} = R(\theta)(m\vec{m}_1 + n\vec{n}_1), \quad a, b, m, n \in \mathbb{Z}
\]

<p align="center">
  <img src="/images/moire.jpeg" width="400"/>
</p>

For a commensurate moiré superlattice to form, there must exist integers \( a, b, m, n \) such that:

\[
|a\vec{a}_1 + b\vec{b}_1| = |m\vec{m}_1 + n\vec{n}_1|
\]

This leads to a condition based on vector magnitudes and orientations, expressed as:

\[
\vec{a}_1 \cdot \vec{b}_1 = \vec{m}_1 \cdot \vec{n}_1 = \cos{\beta} = \frac{m^2 + n^2 - a^2 - b^2}{2(ab - mn)}
\tag{1}
\label{Diiii}
\]

Given a pair \( (a, b) \in \mathbb{Z} \), the objective is to find integer values \( (m, n) \) that satisfy Eqn. (1). Once such a solution is found, the corresponding twist angle \( \theta \) can be computed.

<br><br>

The procedure for computing \( \theta \) proceeds as follows. The length of the moiré lattice vector \( \vec{r} \) connecting equivalent lattice points is:

\[
r = |\vec{r}| = \sqrt{a^2 + b^2 + 2ab\cos\beta}
\]

Using the Law of Sines, the intermediate angles \( \alpha \) and \( \gamma \) in the triangle are calculated as:

\[
\frac{b}{\sin\alpha} = \frac{r}{\sin(180^\circ - \beta)} \quad \Rightarrow \quad 
\alpha = \sin^{-1}\left(\frac{b \sin\beta}{r}\right)
\tag{2}
\label{one}
\]

\[
\frac{n}{\sin\gamma} = \frac{r}{\sin(180^\circ - \beta)} \quad \Rightarrow \quad 
\gamma = \sin^{-1}\left(\frac{n \sin\beta}{r}\right)
\tag{3}
\label{two}
\]

Finally, the twist angle \( \theta \) between the two layers is given by the difference between these two angles:

\[
\theta = \alpha - \gamma = \sin^{-1}\left(\frac{b \sin\beta}{r}\right) - \sin^{-1}\left(\frac{n \sin\beta}{r}\right)
\tag{4}
\]

This classical approach, based on solving Diophantine conditions, provides a rigorous framework to determine commensurate twist angles where the lattice vectors of the two layers form periodic overlaps.

</div>



<!--A popular approach is to formulate the problem algebraically by equating integer linear combinations of lattice vectors under rotation. That is, you search for integers \( m_1, m_2, n_1, n_2 \) such that:

$$
m_1\mathbf{a}_1 + m_2\mathbf{a}_2 = R(\theta)(n_1\mathbf{b}_1 + n_2\mathbf{b}_2)
$$

This is essentially a **Diophantine equation** in four variables, constrained by a rotation angle \( \theta \).-->

<div style="text-align: justify;">
<div style="height: 10px;"></div>
<h3>Time Complexity of Diophantine Approach</h3>

Let the radius cutoff be r. Then, each lattice contains approximately \( n = O(r^2) \) points. As a result, the Diophantine approach requires checking \( O(n^2) \) combinations. Thus, this becomes increasingly slow as the radius grows. However, our motivation for seeking a simpler and faster method wasn't primarily performance — we just found the Diophantine approach inelegant and unsatisfying.

</div>
<div style="height: 10px;"></div>
## Observations About Lattice Structure

We noticed some very helpful properties in regularly spaced lattices (like triangular or square lattices).

{Insert image showing how points are grouped into concentric distance shells.}

If we arrange all the points in the lattice in ascending order of their distances from the origin, we can clearly see discrete levels forming, sort of like a step function. Each level — meaning, each radius from the origin where lattice points lie — contains **symmetrically placed points**. For example, in a triangular lattice, each such level contains a multiple of 6 points due to 6-fold symmetry.

{Insert plot of number of points per radius level showing multiples of 6.}

This means that if two lattices share the same level (same distance from origin), then **we only need to align one pair of points on that level**, and the rest 6 will align automatically due to symmetry.

## Our Method

We approach the angle-finding task using the geometry and symmetry of lattices, avoiding algebraic equations altogether.

### Algorithm Overview

Let \( A_r \) and \( B_r \) be the sets of lattice points (from lattice A and B respectively) within radius \( r \) from the origin.

1. **Group points by radius**:  
   For each point in \( A_r \) and \( B_r \), compute its distance \( d = \|\mathbf{p}\| \) from the origin.  
   Group points that lie at the same radius (up to numerical tolerance) into *levels*.

2. **Identify shared levels**:  
   Let \( D = \{d \mid d \text{ occurs in both } A_r \text{ and } B_r \} \setminus \{0\} \).  
   These are the radii at which both lattices have points.

3. **Filter by angular sector**:  
   For each \( d \in D \), consider only those points \( \mathbf{p} \in A_r \) and \( \mathbf{q} \in B_r \) on level \( d \) such that  
   $$
   0 < \angle(\mathbf{p}) < \theta_\text{max}, \quad
   0 < \angle(\mathbf{q}) < \theta_\text{max}
   $$  
   where \( \theta_\text{max} \) is the lattice's symmetry sector (e.g., \( 60^\circ \) for triangular lattices).

4. **Compute angle differences**:  
   For all pairs \( (\mathbf{p}, \mathbf{q}) \) from the filtered sets on each shared level \( d \):
   - Let \( \theta_1 = \angle(\mathbf{p}) \), \( \theta_2 = \angle(\mathbf{q}) \)
   - If \( \theta_2 > \theta_1 \) and \( \theta = \theta_2 - \theta_1 > \text{tolerance} \), then:
     - Store angle \( \theta \) along with the point pair \( (\mathbf{p}, \mathbf{q}) \), **only if** no existing pair has a smaller \( \theta_1 \)

{Insert illustration showing concentric circles with marked angular slices and example point pairs.}

This procedure ensures we collect unique, minimal-angle configurations that could align under rotation, constrained to the symmetry of the lattice.

## Time Complexity

Let’s denote:

- \( n = O(r^2) \): number of points in the circular cutoff

### Breakdown:

- **Sorting distances**: \( O(n \log n) \)
- **Finding common radii**: \( O(n) \)
- **Pairwise angle checks**: Very few (multiples of 6) per level → negligible

### Final Time:

$$
O(n \log n)
$$

Compared to the \( O(n^2) \) of Diophantine, this is clearly more efficient — and arguably more intuitive.

## Summary

We avoided solving Diophantine equations by leaning on geometry and symmetry:

- Points are grouped by radius
- Only overlapping radii are considered
- Pairwise angle differences yield the commensurate angles

This makes **MoirePy**’s angle detection both **fast and visual**, and opens up room for further optimizations or generalizations.
