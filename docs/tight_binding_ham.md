# Tight Binding Moiré Hamiltonian Construction

The **tight-binding Hamiltonian** [1] is a widely used model in solid-state physics and quantum chemistry to describe the electronic structure of solids — especially in crystals and layered materials. In this model, electrons are considered *localized* around atomic sites but can *hop* to neighbouring atoms. To describe such a Hamiltonian for a **Moiré system**, we use the second quantized [2] form:

$$
H = \sum_{\alpha, \beta;\, r,r' \in L} t^1_{rr', \alpha\beta}c^{\dagger}_{r,\beta}c_{r',\alpha} 
+ \sum_{\alpha, \beta;\, r,r' \in U} t^2_{rr', \alpha\beta}d^{\dagger}_{r,\beta}d_{r',\alpha} 
+ \sum_{\alpha, \beta;\, r,r'} t^{\perp}_{rr', \alpha\beta}c^{\dagger}_{r,\beta}d_{r',\alpha} + \text{h.c.}
$$

Here, $c^{\dagger}_{r,\beta}$ and $c_{r',\alpha}$ denote the electron creation and annihilation operators at lattice sites $r$ and $r'$ in the **lower layer** ($L$), associated with orbitals $\beta$ and $\alpha$, respectively. Likewise, $d^{\dagger}_{r,\beta}$ and $d_{r',\alpha}$ are the corresponding operators in the **upper layer** ($U$).

The terms $t^1_{rr', \alpha\beta}$ and $t^2_{rr', \alpha\beta}$ represent the *intralayer hopping amplitudes*, describing electron tunneling from orbital $\alpha$ at site $r'$ to orbital $\beta$ at site $r$ within the lower and upper layers, respectively. In the special case where $r = r'$ and $\alpha = \beta$, these terms correspond to the *on-site potential*—the energy of an electron residing in a particular orbital. The *interlayer coupling* is described by $t^{\perp}_{rr', \alpha\beta}$, which governs the hopping of an electron from orbital $\alpha$ at site $r'$ in the upper layer to orbital $\beta$ at site $r$ in the lower layer.

For simplicity, consider only **nearest-neighbour hopping** with a single orbital per site (**MoirePy** can handle any arbitrary number of orbital systems). In such cases, the orbital indices \( \alpha \) and \( \beta \) can be omitted to simplify the notation. We can define the basis as:

$$
\Psi^{\dagger} = (c^{\dagger}_{1}, c^{\dagger}_{2}, \dots, c^{\dagger}_{n}, d^{\dagger}_{1}, d^{\dagger}_{2}, \dots, d^{\dagger}_{n})
$$

Here, \( c^{\dagger}_{i} \) (\( d^{\dagger}_{i} \)) is the creation operator at site \( i \) in the lower (upper) layer.

Then, the Hamiltonian takes the block matrix form:

$$
H = \Psi^{\dagger}
\begin{pmatrix}
h_{LL} & h_{LU} \\
h_{UL} & h_{UU}
\end{pmatrix}
\Psi
$$

Here, \( h_{LL} \) and \( h_{UU} \) are the *first-quantized* \( n \times n \) Hamiltonians of the lower and upper layers, respectively. The blocks \( h_{LU} \) and \( h_{UL} \) represent interlayer couplings.

## References

1. Neil Ashcroft, and David Mermin, Solid State Physics, Saunders College Publishing, 1976.

2. H. Bruus and K. Flensberg, Many-Body Quantum Theory in Condensed Matter Physics: An Introduction, OUP publication, 2004.