# Theoretical Background

The discovery that twisting two atomically-thin layers, like graphene, by a small angle can induce superconductivity and other correlated phases has transformed condensed matter physics. At the heart of this transformation lies the **moiré pattern**: a long-wavelength lattice emerging from the interference of two periodic structures. These patterns reshape the electronic landscape, giving rise to flat bands, topological states, and strong interaction effects absent in the individual layers.

Simulating such systems demands a careful balance between **atomistic resolution** and **large-scale periodicity**. Full ab-initio methods become computationally infeasible at moiré scales, while naive tight-binding models often miss crucial geometrical and symmetry constraints. **MoirePy** is designed to bridge this gap, converting clean geometric input into efficient, flexible, and physically grounded tight-binding models for commensurate moiré lattices.

## Purpose of This Section

This section documents the theoretical models and computational strategies implemented in MoirePy. From the construction of **tight-binding Hamiltonians** and their **real- and $k$-space formulations**, to **geometry-based neighbor search algorithms** and **quantized twist-angle calculations**, we outline the core principles underlying the package.

We believe in **transparent scientific software**, and this section is written not just for end users, but also for researchers and developers who want to understand, modify, or build upon MoirePy. Whether you're reviewing our assumptions, adapting our algorithms, or simply learning the methods behind moiré modeling, these notes aim to make the theoretical foundations as clear and open as possible.


