# Angle Calculation Process

In this section, we describe how **MoirePy** determines the *commensurate rotation angles* between two lattices clipped to a bounded region. Our goal is to compute the angles where two lattices can overlap nicely, giving rise to observable moiré patterns.

## Problem Definition

Let lattice \( A \) and lattice \( B \) be two periodic 2D point sets defined by their lattice vectors.

We are interested in the following:

> Given a rotation \( \theta \), does there exist a point \( \mathbf{p} \in A \) and a point \( \mathbf{q} \in B \) such that  
> 
> $$
\mathbf{p} = R(\theta)\mathbf{q}
$$

Where \( R(\theta) \) is the standard 2D rotation matrix:

$$
R(\theta) = 
\begin{bmatrix}
\cos\theta & -\sin\theta \\
\sin\theta & \cos\theta
\end{bmatrix}
$$

To make the computation finite, we consider only the lattice points within a circular cutoff of radius \( r \). Let’s denote the sets:

- \( A_r = A \cap \text{circle}(r) \)
- \( B_r = B \cap \text{circle}(r) \)

We now want to find angles \( \theta \) such that:

$$
\exists \mathbf{p} \in A_r, \mathbf{q} \in B_r: \quad \mathbf{p} = R(\theta)\mathbf{q}
$$

## Traditional Diophantine Equation Approach

A popular approach is to formulate the problem algebraically by equating integer linear combinations of lattice vectors under rotation. That is, you search for integers \( m_1, m_2, n_1, n_2 \) such that:

$$
m_1\mathbf{a}_1 + m_2\mathbf{a}_2 = R(\theta)(n_1\mathbf{b}_1 + n_2\mathbf{b}_2)
$$

This is essentially a **Diophantine equation** in four variables, constrained by a rotation angle \( \theta \).

### Time Complexity

Let the radius cutoff be \( r \). Each lattice will have approximately \(n = O(r^2) \) points.

The Diophantine approach needs to check \( O(n^2) \) combinations, since both lattices have n points.

This gets very slow as radius increases. We wanted a simpler and faster method — not because this is slow, but because we just didn’t like it.

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
