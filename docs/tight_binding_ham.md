## Tight Binding Moiré Hamiltonian Construction

The **tight-binding Hamiltonian** is a widely used model in solid-state physics and quantum chemistry to describe the electronic structure of solids — especially in crystals and layered materials. In this model, electrons are considered *localized* around atomic sites but can *hop* to neighboring atoms. To describe such a Hamiltonian for a **Moiré system**, we use the second quantized form:

$$
H = \sum_{\alpha, \beta;\, r,r' \in L} t^1_{rr', \alpha\beta}c^{\dagger}_{r,\beta}c_{r',\alpha} 
+ \sum_{\alpha, \beta;\, r,r' \in U} t^2_{rr', \alpha\beta}d^{\dagger}_{r,\beta}d_{r',\alpha} 
+ \sum_{\alpha, \beta;\, r,r'} t^{\perp}_{rr', \alpha\beta}c^{\dagger}_{r,\beta}d_{r',\alpha} + \text{h.c.}
$$

Here, 
* $c^{\dagger}_{r,\beta}$, ($c_{r',\alpha}$): creation (annihilation) operators at sites $r$, $r'$ for orbitals $\beta$, $\alpha$ in the **lower layer** $L$.
* $d^{\dagger}_{r,\beta}$, $d_{r',\alpha}$: corresponding operators in the **upper layer** $U$.
* $t^1_{rr', \alpha\beta}$, $t^2_{rr', \alpha\beta}$: hopping amplitudes within the lower and upper layers, respectively.
* $t^{\perp}_{rr', \alpha\beta}$: interlayer hopping amplitude from orbital $\alpha$ at $r'$ in the upper layer to orbital $\beta$ at $r$ in the lower layer.
* When $r = r'$ and $\alpha = \beta$, the hopping term reduces to the **on-site potential**.

---

### Simplified Case: Nearest-Neighbor Hopping

For many purposes, we can simplify the model by assuming:

* Only **nearest-neighbor hopping**, and
* A **single orbital per site**.

In this case, we can omit the orbital indices $\alpha$, $\beta$, and write the Hamiltonian in matrix form.

Let the basis be:

$$
\Psi^{\dagger} = \left(c^{\dagger}_{1}, c^{\dagger}_{2}, \dots, c^{\dagger}_{n}, d^{\dagger}_{1}, d^{\dagger}_{2}, \dots, d^{\dagger}_{n} \right)
$$

Then the Hamiltonian becomes:

$$
H = \Psi^{\dagger}
\begin{pmatrix}
h_{LL} & h_{LU} \\
h_{UL} & h_{UU}
\end{pmatrix}
\Psi
$$

Where:

* $h_{LL}$ and $h_{UU}$: first-quantized $n \times n$ Hamiltonians for the lower and upper layers.
* $h_{LU}$ and $h_{UL}$: interlayer coupling blocks.

In **Chapter 2**, we will explore how to construct such Hamiltonians using the `MoirePy` library with minimal effort.

---

## K-space Hamiltonian

If the system exhibits **translational invariance**, then **momentum** (or crystal momentum) becomes a good quantum number. In that case, we can **diagonalize** the Hamiltonian in momentum space, making it easier to:

* Obtain the **band structure**, and
* Visualize the energy spectrum as a function of momentum.

In the next section, we will explain how to obtain the **K-space Hamiltonian**, i.e., the Hamiltonian expressed in the **momentum eigenbasis**.

---

Let me know if you'd like this exported to a `.md` file or adapted for a specific renderer like Jupyter or MkDocs.




### Rerefences

