# MoirePy: Twist It, Solve It, Own It!

**MoirePy** is a FOSS Python library for simulating **moire lattices** with a clean, highly Pythonic API. The performance-critical backend is written in **Rust**, which keeps large lattice generation and Hamiltonian assembly fast while preserving a Python-first workflow.

It is designed to be easy to use while remaining fully flexible. You can define custom lattice geometries, program arbitrary hopping functions, and build tight-binding Hamiltonians in both **real space** and **k-space**, with support for **open (OBC)** and **periodic (PBC)** boundary conditions. Generate your lattice and Hamiltonian in MoirePy, then export to the tools you already use (like **Kwant**, **NumPy**, **SciPy**, and related ecosystems) and keep your workflow.

* **Documentation:** [https://jabed-umar.github.io/MoirePy/](https://jabed-umar.github.io/MoirePy/)<br>
* **Github Repository:** [https://github.com/jabed-umar/MoirePy](https://github.com/jabed-umar/MoirePy)<br>
* **PyPI page:** [https://pypi.org/project/moirepy/](https://pypi.org/project/moirepy/)

<p align="center">
  <!-- Replace with GIF later -->
  <img src="https://jabed-umar.github.io/MoirePy/images/cover_image_gemini.webp" width="600">
</p>

---

## Why MoirePy

* **Plug and play**: Build in MoirePy, continue in your current stack: Kwant, NumPy, SciPy etc.
* **Fast by default**: Rust backend + optimized core algorithms.
* **Flexible models**: Custom geometry, orbitals, and hoppings.
* **Multiple modes**: Real-space/k-space and OBC/PBC supported.
* **Easy onboarding**: Pythonic API and web tools for quick tinkering.

---

## Philosophy

* **Not a workflow replacement**: Use MoirePy where it shines, then export to *Kwant*, *NumPy*, *SciPy*, and more.
* **Do one thing and do it well**: Focus on robust and fast moire lattice simulation and Hamiltonian generation.
* **Power to you**: If it is mathematically possible, you can build it. No questions asked.
* **Learn by doing**: Explore fast, and reproduce well-known paper systems.

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

## Benchmark: Why This Is Fast

Basic benchmark snapshots using Twisted Bilayer Graphene (TBG) as a test case, benchmarked on a laptop with a 12th Gen Intel Core i5 (performance cores).

### Experiment 1: Full Pipeline

![Experiment 1: Full pipeline benchmark](https://jabed-umar.github.io/MoirePy/benchmark/experiment1_total_pipeline.webp)

This includes lattice generation plus Hamiltonian assembly.
Even at around **100k lattice points**, the full pipeline is about **400 ms**. For context, configurations near the magic angle ($\sim 1.1^\circ$) typically have just around **5-6k lattice points**.

### Experiment 2: Subsequent Hamiltonian Generation

![Experiment 2: Hamiltonian-only benchmark](https://jabed-umar.github.io/MoirePy/benchmark/experiment2_hamiltonian_only.webp)

After setup, Hamiltonian generation is much faster.
At around **100k lattice points**, it is about **30 ms**.
This is the speed profile we optimize for: heavy setup once, then fast repeated builds.
Memory usage in these runs is negligible relative to typical laptop RAM.

For more comprehensive performance benchmarks, visit this: [BENCHMARK](https://jabed-umar.github.io/MoirePy/benchmarks/)

---

## Contributing

Contributions are welcome.

Please read the full contribution guide at [contributing.md](contributing.md).

Highlights:

* Join our [Discord server](https://discord.gg/EzRzMXgzfe) to coordinate and discuss ideas.
* Choose your path: **Docs Contributor** (examples, blogs, docs) or **Developer Contributor** (Rust core, Python API, tests/CI).
* Start from issues and submit focused pull requests.
* Follow project philosophy: keep implementations simple, benchmark in practice, and optimize for user experience.

---

## Citation

If you use this work in research:

```bibtex
@software{MoirePy2025,
  author = {Mukhopadhyay, Aritra and Umar, Jabed},
  license = {MIT},
  title = {{MoirePy: Python package for efficient atomistic simulation of moiré lattices}},
  url = {https://github.com/jabed-umar/MoirePy}
}
```
