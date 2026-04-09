---
title: 'MoirePy: A Python-Rust Framework for Efficient Simulation of Moiré Lattices'
tags:
  - Python
  - Rust
  - Moiré lattices
  - Condensed matter physics
  - Tight-binding
  - Twisted bilayer graphene
authors:
  - name: Aritra Mukhopadhyay
    orcid: 0009-0005-7960-373X
    affiliation: "1, 3, 4"
    corresponding: true
    email: aritra.mukhopadhyay@niser.ac.in
  - name: Jabed Umar
    orcid: 0009-0000-0566-4373
    affiliation: "2, 3, 4"
    corresponding: true
    email: aritra.mukhopadhyay@niser.ac.in
affiliations:
  - name: "Independent Researcher, India"
    index: 1
  - name: "Technischen Universität Wien (TU Wien), Vienna, Austria"
    index: 2
  - name: "National Institute of Science Education and Research (NISER), Bhubaneswar, India"
    index: 3
  - name: "Homi Bhabha National Institute (HBNI), Mumbai, India"
    index: 4
date: 25 March 2026
bibliography: paper.bib
---

# Summary

MoirePy is a Python library with a Rust backend for high-performance simulation of moiré superlattices. It provides a complete pipeline for constructing commensurate twisted bilayer systems, generating atomic coordinates, and assembling tight-binding Hamiltonians. The library supports arbitrary lattice geometries, user-defined hopping functions, and both open and periodic boundary conditions. Hamiltonians are returned as SciPy sparse matrices [@virtanen2020scipy], making MoirePy interoperable with NumPy [@harris2020array], Kwant [@groth2014kwant], and the broader scientific Python ecosystem.

# Statement of Need

The discovery of unconventional superconductivity in magic-angle twisted bilayer graphene [@cao2018unconventional; @cao2018correlated] has driven intense research into moiré systems. A central challenge is computational cost: near the magic angle ($\sim 1.1°$), commensurate unit cells contain thousands of atoms, and exploring the twist-angle parameter space requires constructing and solving many such systems. While tools like Kwant [@groth2014kwant] excel at quantum transport calculations, there is currently no convenient and fast way to go from a twist angle to a ready-to-use moiré Hamiltonian. Researchers must handle lattice construction, neighbour searching, and Hamiltonian assembly largely by hand, often stitching together ad-hoc scripts that are slow for large systems. MoirePy was built to fill this gap: a single, high-performance library that automates the entire pipeline and exports standard sparse matrices for downstream analysis.

# State of the Field

Existing tools like Twister [@twister2022], MoireStudio [@moirestudio2026], and PyMoire [@pymoire] often attempt to cover the entire simulation pipeline, leading to "half-done" features and rigid workflows. Instead of creating another incomplete, all-in-one suite, MoirePy adopts the Unix philosophy: doing one thing well. We focus exclusively on the efficient generation of lattice geometries and Hamiltonians. By providing a high-performance, Rust-backed engine that exports standard sparse matrices to specialized analysis software like Kwant [@groth2014kwant], MoirePy serves as a modular building block optimized for speed and interoperability rather than a monolithic framework.

# Software Design

MoirePy is designed around two principles: make the common tasks convenient, and make the heavy computations fast.

The first bottleneck in any moiré simulation is selecting a valid commensurate twist angle. MoirePy includes a geometric shell-based algorithm that enumerates all commensurate angles for a given lattice type. To make this step even more accessible, the project provides an interactive, web-based Angle Value Calculator where users can explore the angle space and pick configurations without writing any code.

The second design goal is a clean separation between geometry and physics. The lattice and its inter-site connections are constructed once; the user then generates Hamiltonians from that geometry as many times as needed, each time with different hopping parameters. Building the connections requires a spatial neighbour search, which scales as $O(N \log N)$ for $N$ lattice sites using a KD-tree [@kiddo]. Once the connections are established, each subsequent Hamiltonian assembly runs in $O(N)$ time. Hopping amplitudes can be simple constants or arbitrary callable functions, making it straightforward to implement distance-dependent models such as Slater-Koster parametrisations [@slater1954simplified] or Gaussian decay without touching the geometry code.

Beyond algorithmic choices, all performance-critical operations (angle enumeration, neighbour search, connectivity construction) are implemented in Rust and exposed to Python via PyO3. This combination of efficient algorithms and a compiled backend is what produces the benchmark results below.

# Benchmarks

We benchmarked MoirePy on twisted bilayer graphene systems of increasing size on a laptop with a Intel Core i5 12th gen processor. All the experiemeents have been repeated 100 times and an average of the results is presented.

![Full pipeline benchmark for twisted bilayer graphene on a hexagonal lattice. Even at 100k lattice sites the total time remains under 400 ms.\label{fig:exp1}](benchmark/experiment1_total_pipeline.pdf){ width=80% }

\autoref{fig:exp1} shows the full pipeline (lattice points and connection generation plus one time Hamiltonian assembly). A system with approximately 100k lattice sites completes in less than 400 ms. Configurations near the magic angle typically involve only 5-6k sites and finish in under 15 ms.

![Subsequent Hamiltonian generation benchmark. After initial lattice setup, regenerating the Hamiltonian at 100k sites takes only 12 ms.\label{fig:exp2}](benchmark/experiment2_hamiltonian_only.pdf){ width=80% }

\autoref{fig:exp2} isolates subsequent Hamiltonian generation after the geometry has been constructed. This is the relevant metric for parameter sweeps where hopping values change but geometry stays fixed. At 100k sites, regeneration takes approximately 12 ms. Where as general use cases involving a few thousand sites is completed in under 0.25 ms.

These results demonstrate that MoirePy handles systems well beyond magic-angle scales while keeping iteration times short enough for interactive exploration on commodity hardware.

# Research Impact

MoirePy accelerates moiré quantum matter research by bridging geometric construction with high-performance solvers. Its impact is centered on three areas:

* **High-Throughput Exploration:** Electronic properties in twisted bilayer graphene (TBG) are hypersensitive to the twist angle $\theta$. Identifying flat-band regimes requires scanning dense parameter spaces. MoirePy reduces Hamiltonian assembly from minutes to milliseconds, enabling exhaustive scans on standard workstations that previously required HPC clusters.
* **Versatility in Boundary Physics:** A unified API supports both Periodic (PBC) and Open Boundary Conditions (OBC). This facilitates rapid switching between reciprocal-space calculations—essential for band structures and Chern numbers—and real-space investigations of structural relaxation, impurities, and finite-size effects in moiré "flakes."
* **Ecosystem Interoperability:** By outputting standard SciPy sparse matrices, MoirePy acts as a modular pre-processor for the scientific Python stack. Researchers can immediately leverage established solvers like **Kwant** for transport [@groth2014kwant], **PySCF** for mean-field theory, or **PETSc** for large-scale eigenproblems, avoiding the "walled garden" of monolithic tools.

# Acknowledgements
The authors are grateful to Sarbajit Mazumdar (Julius-Maximilians-Universität Würzburg) for insightful discussions on the physics of moiré systems and for validating the library through the replication of seminal literature. We also acknowledge the National Institute of Science Education and Research (NISER), Bhubaneswar, for providing the environment where this collaboration first began.

# AI Usage Disclosure

AI tools (like GitHub Copilot) assisted in refactoring code and polishing prose. All technical content, core design decisions, and benchmark methodology were produced by the human authors. The final manuscript was reviewed and validated by the authors.

# References

