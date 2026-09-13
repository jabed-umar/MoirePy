# MoirePy API Design Brainstorming

This document logs the independent API proposals from the blind Python API Architect for MoirePy.

## Phase 1: Blind Proposal (Python API Architect)

**Philosophy**: SciPy/Kwant modularity. Separating Geometry from Physics.
**Draft**:
```python
# 1. Geometry (Support for Heterobilayers natively)
bottom = mp.HexagonalLattice(a=2.46, name="bottom")
top = mp.HexagonalLattice(a=2.46, name="top")
twist = mp.CommensurateTwist(ll1=31, ll2=32, ul1=32, ul2=31)
bilayer = mp.MoireBilayer(bottom, top, twist_params=twist)

# 2. Physics (Tight Binding Builder)
model = mp.TightBinding(bilayer)
model.add_intralayer_hopping("bottom", max_dist=1.5, t=-2.7)
model.set_interlayer_hopping("bottom", "top", func=lambda d: 0.3 * np.exp(-d))

# 3. Solvers
H_real = model.hamiltonian(space="real", boundary="pbc")
H_k = model.hamiltonian(space="momentum")
```

## Phase 2: Verdict vs Current MoirePy API

**The Harsh Truth**: The blind agent's API is significantly better than the current MoirePy API.

**Why the current MoirePy API is lacking:**
1. **Monolithic Constructor**: Currently, `BilayerMoireLattice` takes `latticetype=HexagonalLayer` and applies it to *both* layers. You cannot easily do Heterobilayers (e.g., Graphene on hBN). The blind agent fixes this by instantiating `bottom` and `top` separately.
2. **Hardcoded Physics**: Currently, `lattice.generate_hamiltonian(tll=1, tuu=1, tlu=1, tul=1)` hardcodes exactly 4 nearest-neighbor hopping parameters. This is extremely rigid. What if a user wants next-nearest neighbor? What if they want a continuous distance-decay function for interlayer hopping? The blind agent solves this by using a `TightBinding` builder class that accepts distance cutoffs and callable functions.
