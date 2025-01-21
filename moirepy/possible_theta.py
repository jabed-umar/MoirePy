import numpy as np

from moirepy.layers import (
    Layer,
    HexagonalLayer,
    SquareLayer,
    RhombusLayer,
    TriangularLayer,
    KagomeLayer,
)

from tqdm import tqdm, trange


def check_abmn(a, b, m, n, a1, b1, a2, b2, tol):
    if np.abs(round(n) - n) > tol:
        return None
    n = round(n)
    
    if a == 0 or b == 0 or a == m:
        return None
    
    # find the angle between (a*a1 + b*b1) and (m*a2 + n*b2)
    vec1 = a*a1 + b*b1
    vec2 = m*a2 + n*b2
    
    cos = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    
    return a, b, m, n, float(np.rad2deg(np.arccos(cos)))

def find_values(layer1: HexagonalLayer, layer2: HexagonalLayer, limit: int = 20, tol: float = 1e-4) -> None:
    layer1, layer2 = layer1(), layer2()
    
    a1 = layer1.lv1
    b1 = layer1.lv2
    
    a2 = layer2.lv1
    b2 = layer2.lv2
    
    # cos angle between a and b
    cos1 = np.dot(a1, b1) / (np.linalg.norm(a1) * np.linalg.norm(b1))
    cos2 = np.dot(a2, b2) / (np.linalg.norm(a2) * np.linalg.norm(b2))
    abs_n =  np.linalg.norm(b2)
    print(cos1, cos2, abs_n)
    
    found = []
    for a in trange(1, limit+1):
        for b in range(1, limit+1):
            for m in range(1, limit+1):
                
                aa = a * np.linalg.norm(a1)
                bb = b * np.linalg.norm(b1)
                mm = m * np.linalg.norm(a2)

                # solving using sridhar acharya
                sa_b = 2 * mm * cos2
                sa_c = mm*mm - (aa*aa + bb*bb + 2*aa*bb*cos1)
                # print("", sa_b, sa_c, aa, bb, mm, abs_n, sep="    ")

                if sa_b*sa_b - 4*sa_c < 0: continue  # no real solution
                
                
                n = (-sa_b + np.sqrt(sa_b*sa_b - 4*sa_c)) / (2*abs_n)
                res = check_abmn(a, b, m, n, a1, b1, a2, b2, tol)
                if res:
                    found.append(res)

                n = (-sa_b - np.sqrt(sa_b*sa_b - 4*sa_c)) / (2*abs_n)
                res = check_abmn(a, b, m, n, a1, b1, a2, b2, tol)
                if res:
                    found.append(res)
    
    # sort the found array according to the last Values (theta)
    found.sort(key=lambda x: x[-1])
    
    return found


if __name__ == "__main__":
    # matches = find_values(SquareLayer, SquareLayer)
    matches = find_values(RhombusLayer, RhombusLayer, limit=20)
    print(len(matches))
    print(*matches, sep="\n")
