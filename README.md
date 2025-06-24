# MoirePy: Twist It, Solve It, Own It!

**MoirePy** is a FOSS Python package for the simulation and analysis of **bilayer moiré lattices** using **tight-binding models**. Built for computational physicists and material scientists, it enables quick and flexible moiré band structure calculations, visualization, and manipulation. Our primary focus is on **commensurate** moiré lattices only.


**Documentation:** [https://jabed-umar.github.io/MoirePy/](https://jabed-umar.github.io/MoirePy/)<br>
**Github Repository:** [https://github.com/jabed-umar/MoirePy](https://github.com/jabed-umar/MoirePy)<br>
**PyPI page:** [https://pypi.org/project/moirepy/](https://pypi.org/project/moirepy/)

## Features

- Fast and efficient simulation of 2D bilayer moiré lattices.
- Efficient $O(\log n)$ time nearest neighbour searches.
- supports **custom lattice definitions** with some basic predefined ones:
    - Triangular
    - Square
    - Hexagonal
    - Kagome
- both **real** and **k-space Hamiltonian** generation for tight-binding models with:
    - Nearest-neighbour coupling
    <!-- - Nth nearest-neighbour coupling -->
    - Arbitrary number of orbitals per site
    - All couplings can be real (default), or complex numbers.
    - All couplings can be functions of position of the point(s) and the point type(s) (for example, different coupling for A-A, A-B, B-B sites for hexagonal lattices)
    - Custom Intra and Interlayer Coupling Design.
- [Web based tool](https://jabed-umar.github.io/MoirePy/theory/avc/) makes it convenient to calculate lattice angles before simulation.
- Extensive Documentation and examples for easy onboarding.
- Compatible with other related libraries like Kwant (so that you can generate moire Hamiltonian and use it with Kwant for further analysis).
- **Freedom to researcher:** We allow you to define your layers and apply whatever couplings you want. If you want the lattice points to have 53 orbitals each—sure, go ahead. As long as you know what you're doing, we won’t stop you. We don't verify whether it's physically possible.

## Upcoming Features

- **Support for higher-dimensional layers**: Extend current 2D-only support to include higher dimensional constituent layers.
- **Multi-layer stacking**: Go beyond bilayers; enable simulation of trilayers and complex heterostructures.
- **Non-equilibrium Green's function support** *(research in progress)*: Develop tools for computing Green’s functions efficiently to study non-equilibrium and quantum transport phenomena.

## Installation

You can install MoirePy from PyPI via pip:

```bash
$ pip install moirepy
```

## Basic Usage

For detailed usage, please refer to our [documentation](https://jabed-umar.github.io/MoirePy/).

```python
>>> import matplotlib.pyplot as plt
>>> from moirepy import BilayerMoireLattice, TriangularLayer
>>> # Define the Moiré lattice with two triangular layers
>>> moire_lattice = BilayerMoireLattice(
>>>     latticetype=TriangularLayer,
>>>     ll1=9, ll2=10,
>>>     ul1=10, ul2=9,
>>>     n1=1, n2=1,  # number of unit cells
>>> )
twist angle = 0.0608 rad (3.4810 deg)
271 points in lower lattice
271 points in upper lattice
>>> ham = moire_lattice.generate_hamiltonian(
>>>     tll=1, tuu=1, tlu=1, tul=1,
>>>     tuself=1, tlself=1,
>>> )
>>> plt.matshow(ham, cmap="gray")
```
![alt text](docs/images/getting_started/moire_lattice_and_hamiltonian/out2.png)

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



## Cite This Work

If you use this software or a modified version in academic or scientific research, please cite:

```BibTeX
@misc{MoirePy2025,
	author = {Aritra Mukhopadhyay, Jabed Umar},
	title = {MoirePy: Python package for efficient atomistic simulation of moiré lattices},
	year = {2025},
	url = {https://jabed-umar.github.io/MoirePy/},
}
```
