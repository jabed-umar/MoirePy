"""MoirePy v2 - A pure Python, N-layer tight-binding Moire stack builder."""

from .utils import (
    get_rotation_matrix,
    are_coeffs_integers,
    LatticeAlreadyFinalisedError,
    StackStateError,
)
from .layers import Layer

__all__ = [
    "get_rotation_matrix",
    "are_coeffs_integers",
    "LatticeAlreadyFinalisedError",
    "StackStateError",
    "Layer",
]
from .presets import (
    SquareLayer,
    TriangularLayer,
    HexagonalLayer,
    KagomeLayer,
    Rhombus60Layer,
)

__all__.extend([
    "SquareLayer",
    "TriangularLayer",
    "HexagonalLayer",
    "KagomeLayer",
    "Rhombus60Layer",
])
from .stack import MoireStack
__all__.append("MoireStack")
