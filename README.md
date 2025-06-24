# MoirePy: Twist It, Solve It, Own It!

The **moiré effect** is a fascinating visual phenomenon seen in *modern art*, *textile patterns*, and even *currency anti-counterfeiting techniques*, where overlapping periodic structures create intricate interference patterns. But beyond its aesthetic appeal, **moiré effect** has taken centre stage in **condensed matter physics**, where geometry meets electronics, and simple twists unlock a playground of quantum phenomena. 

The discovery of emergent phenomena in **moiré materials**—such as *superconductivity*, *correlated insulating states*, and *topological phases of matter* has sparked a revolution in condensed matter physics. These systems, formed by stacking two-dimensional crystals at small twist angles, give rise to rich electronic landscapes governed by long-wavelength moiré patterns. As the field rapidly expands, the need for intuitive, efficient, and high-precision computational tools has never been greater.

**MoirePy** is an open source (free) Python package built to support the needs of *researchers, theorists, computational* and *material scientists* for the numerical calculation of moiré lattices with a focus on commensurate structure. It enables fast, flexible *atomistic simulations of commensurate moiré lattices*, with tools to compute *dispersion relation*, visualize *band structures* and *density of states*, study *non-hermitian systems* and explore *quantum transport* phenomena with precision and control. It can simulate *metal*, *topological insulators*, *quantum hall effect, superconductivity, spintronics*, or any combination. 


**Documentation:** [https://jabed-umar.github.io/MoirePy/](https://jabed-umar.github.io/MoirePy/)
<br>
**Github Repository:** [https://github.com/jabed-umar/MoirePy](https://github.com/jabed-umar/MoirePy)
<br>
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
    - Nth nearest-neighbour coupling
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

<!--## Installation

You can install MoirePy via pip:

```bash
pip install moirepy
```

 ## Basic Usage

For detailed usage, please refer to our [documentation](https://jabed-umar.github.io/MoirePy/).

```python
from moirepy import MoireLattice  # this i think we dont need as we have one dedicated installation page.
#### u might add  a line like cheek here (link inserted) to install the MoirePy
``` -->


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
