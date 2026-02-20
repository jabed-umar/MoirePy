# MoirePy: Twist It, Solve It, Own It!

**MoirePy** is a FOSS Python library for simulating **moire lattices** with a clean, highly Pythonic API.

It is designed to be easy to use while remaining fully flexible. You can define custom lattice geometries, program arbitrary hopping functions, and build tight-binding Hamiltonians in both **real space** and **k-space**, with support for **open (OBC)** and **periodic (PBC)** boundary conditions. Generated Hamiltonians can also be exported to tools like Kwant for further analysis.

* **Documentation:** [https://jabed-umar.github.io/MoirePy/](https://jabed-umar.github.io/MoirePy/)<br>
* **Github Repository:** [https://github.com/jabed-umar/MoirePy](https://github.com/jabed-umar/MoirePy)<br>
* **PyPI page:** [https://pypi.org/project/moirepy/](https://pypi.org/project/moirepy/)

<p align="center">
  <!-- Replace with GIF later -->
  <img src="https://upload.wikimedia.org/wikipedia/commons/7/70/Example.png" width="500">
</p>

---

## Why MoirePy

* **Define anything**: No restrictions on lattice geometry or orbitals
* **Pythonic API**: If you know Python, MoirePy is intuitive
* **Custom hoppings**: Fully programmable intra/inter-layer couplings
* **Real-space + k-space**: Both supported natively
* **OBC + PBC**: Switch boundary conditions easily
* **Kwant-compatible**: Export Hamiltonians for further analysis
* **Fast construction**: Efficient neighbour search using KDTree
* **Reproducibility-first**: Designed to replicate known results and papers

---

## Installation

```bash
pip install moirepy
```

---

## Quick Example: Twisted Bilayer Graphene

```python
import numpy as np
import matplotlib.pyplot as plt
from moirepy import BilayerMoireLattice, HexagonalLayer

# Define a twisted bilayer moiré lattice
lattice = BilayerMoireLattice(
    latticetype=HexagonalLayer,
    # you choose the next 4 values based on the twist angle using this tool:
    # Angle-Value Calculator: https://jabed-umar.github.io/MoirePy/theory/avc/
    ll1=3, ll2=4, ul1=4, ul2=3,
    n1=1, n2=1,
)

# Visualize the lattice
lattice.plot_lattice()
plt.show()

ham = lattice.generate_hamiltonian(
    tll=1, tuu=1,
    tlu=1, tul=1,
    tlself=0, tuself=0
)  # returns scipy sparse matrix

print(ham.shape)
```

---

## Philosophy

MoirePy does not try to enforce what is “physically valid”.

If you want:

* unusual lattices
* non-standard couplings
* high orbital counts

you are free to do so.

> The library gives you control. You decide what makes sense.

---

## Contributing

Contributions are welcome.

* Report bugs or request features via issues
* Submit pull requests for improvements
* Add examples, tutorials, or benchmarks

A detailed contributing guide will be added soon.

---

## Citation

If you use this work in research:

```bibtex
@misc{MoirePy2025,
  author = {Aritra Mukhopadhyay and Jabed Umar},
  title = {MoirePy: Python package for efficient atomistic simulation of moiré lattices},
  year = {2025},
  url = {https://jabed-umar.github.io/MoirePy/}
}
```
