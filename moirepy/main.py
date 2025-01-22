import numpy as np

from layers import (
    Layer,
    HexagonalLayer,
    SquareLayer,
    TriangularLayer,
    RhombusLayer,
    KagomeLayer,
)


class MoireLattice:
    def __init__(self, layer1: Layer, layer2: Layer, a: int, b: int, m: int, n: int) -> None:
        self.layer1 = layer1()
        self.layer2 = layer2()
        
        self.mlv1 = a * self.layer1.lv1 + b * self.layer1.lv2
        self.mlv2 = m * self.layer1.lv1 + n * self.layer1.lv2

        # print(f"{self.mlv1 = }")
        # print(f"{self.mlv2 = }")
        # print(f"{m * self.layer2.lv1 + n * self.layer2.lv2}")
        # print(np.isclose(self.mlv2, m * self.layer2.lv1 + n * self.layer2.lv2))
        
        if not np.isclose(self.mlv2, m * self.layer2.lv1 + n * self.layer2.lv2).all():
            raise ValueError(f"Invalid lattice vectors {a = }, {b = }, {m = }, {n = }. {self.mlv2} != {m * self.layer2.lv1 + n * self.layer2.lv2}")


if __name__ == "__main__":
    layer1 = RhombusLayer
    layer2 = RhombusLayer
    
    lattice = MoireLattice(layer1, layer2, 10, 9, 9, 10)