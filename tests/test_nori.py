import pytest
import numpy as np
import os
import json

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data", "nori")
META_FILE = os.path.join(HERE, "data", "meta.json")

def test_nori(record_mode):
    from moirepy.layers import Layer
    from moirepy.stack import MoireStack

    A_GRAPHENE = 2.46
    D_SPACING = 3.35
    ALPHA = (6.715, 0.762, 0.179, 1.411)
    BETA = (6.811, 0.01, 19.176)
    C_SCALE = 1.44
    T_INTRA = 2.57

    lv1 = np.array([1.0, 0.0]) * A_GRAPHENE
    lv2 = np.array([0.5, np.sqrt(3.0) / 2.0]) * A_GRAPHENE

    basis_points = [
        (0.0, 0.0, "A"),
        (A_GRAPHENE, A_GRAPHENE / np.sqrt(3.0), "B"),
    ]

    neigh_deltas = np.array([
        [0.0, 1.0 / np.sqrt(3.0), 0.0],
        [-0.5, -1.0 / (2.0 * np.sqrt(3.0)), 0.0],
        [0.5, -1.0 / (2.0 * np.sqrt(3.0)), 0.0],
    ]) * A_GRAPHENE

    neighbours = {"A": neigh_deltas, "B": -neigh_deltas}

    l1 = Layer(lv1, lv2, basis_points, neighbours, pbc=True)
    l2 = Layer(lv1, lv2, basis_points, neighbours, pbc=True)

    stack = MoireStack(l1)
    stack.add_layer(l2, 1, 2, 2, 1, translate=(-A_GRAPHENE, -A_GRAPHENE / np.sqrt(3.0), D_SPACING))
    stack.generate_points(n1=1, n2=1)

    r_3d = np.sqrt(D_SPACING**2 + (A_GRAPHENE / (2.0 * np.sqrt(3.0)) + 1e-4)**2)
    stack.generate_connections(np.array([[0, r_3d], [r_3d, 0]]))

    def screening_S(pos_i, pos_j, all_pos):
        rij = np.linalg.norm(pos_i - pos_j)
        ril = np.linalg.norm(all_pos - pos_i, axis=1)
        rjl = np.linalg.norm(all_pos - pos_j, axis=1)
        mask = (ril > 1e-6) & (rjl > 1e-6)
        x = ((ril[mask] + rjl[mask]) / rij).sum()
        xi = BETA[0] * np.exp(-BETA[1] * x**BETA[2])
        return np.tanh(xi)

    def v_sigma(rij):
        r_scaled = C_SCALE * rij
        return ALPHA[0] * r_scaled ** (-ALPHA[1]) * np.exp(-ALPHA[2] * r_scaled ** ALPHA[3])

    def t_perp_func(pos_u, pos_l, all_pos, d):
        rij = np.linalg.norm(pos_u - pos_l)
        inplane = np.linalg.norm(pos_u[:2] - pos_l[:2])
        cos_a = d / np.sqrt(d**2 + inplane**2)
        s_val = screening_S(pos_u, pos_l, all_pos)
        return (cos_a**2) * v_sigma(rij) * (1.0 - s_val)

    def tul(coo_i, coo_j, R_shifts, labels_i, labels_j, stack):
        all_pos = np.vstack([stack.layers[0]._bigger_points, stack.layers[1]._bigger_points])
        vals = []
        for i_idx, j_idx, R in zip(coo_i, coo_j, R_shifts):
            pu = stack.layers[1]._points[i_idx]
            pl = stack.layers[0]._points[j_idx] + R
            vals.append(t_perp_func(pu, pl, all_pos, D_SPACING))
        return np.array(vals)
        
    def t_intra(coo_i, coo_j, R_shifts, labels_i, labels_j, stack):
        return -T_INTRA * np.ones(len(coo_i))

    stack.generate_hamiltonian(
        t_matrix={(1, 0): tul},
        t_self={0: t_intra, 1: t_intra}
    )
    
    mat = np.column_stack([stack.mlv1, stack.mlv2])
    inv_t = 2.0 * np.pi * np.linalg.inv(mat).T
    g1, g2 = inv_t[:, 0], inv_t[:, 1]

    gamma_point = np.zeros(2)
    k1_point = (g1 + 2.0 * g2) / 3.0
    k2_point = (2.0 * g1 + g2) / 3.0

    path_pts = [gamma_point, k1_point, k2_point, gamma_point]
    nseg = 100
    k_list = []
    
    for i in range(len(path_pts) - 1):
        is_last = i == len(path_pts) - 2
        n_pts = nseg + 1 if is_last else nseg
        tvals = np.linspace(0.0, 1.0, n_pts, endpoint=is_last)
        segment = (1.0 - tvals)[:, None] * path_pts[i] + tvals[:, None] * path_pts[i + 1]
        k_list.extend(segment)
        
    k_points = np.array(k_list)
    
    ham_real_space = stack.H.toarray()
    
    all_energies = []
    for kvec in k_points:
        phase = stack.get_phase(kvec)
        ham_k_coo = stack.H.copy()
        ham_k_coo.data = ham_k_coo.data * phase
        ham_k = ham_k_coo.toarray()
        eigvals = np.linalg.eigvalsh(ham_k)
        all_energies.append(eigvals)

    all_energies = np.array(all_energies)
    mid_idx = all_energies.shape[1] // 2
    lower_band = all_energies[:, mid_idx - 1]
    upper_band = all_energies[:, mid_idx]

    v_top = np.max(lower_band)
    c_bottom = np.min(upper_band)
    mu_shift = (v_top + c_bottom) / 2.0

    gap = c_bottom - v_top
    all_energies_shifted = all_energies - mu_shift
    
    current_energies = np.array(all_energies_shifted.tolist(), dtype=np.float16)
    anchor_row = current_energies[75]
    gap_scaled = round(gap * 1e3, 2)
    mu_scaled = round(mu_shift * 1e6, 2)

    if record_mode:
        os.makedirs(DATA_DIR, exist_ok=True)
        anchor_row.tofile(os.path.join(DATA_DIR, "row75.bin"))
        meta = {}
        if os.path.exists(META_FILE):
            with open(META_FILE, "r") as f:
                meta = json.load(f)
        meta["nori"] = {
            "gap_scaled": gap_scaled,
            "mu_scaled": mu_scaled
        }
        with open(META_FILE, "w") as f:
            json.dump(meta, f)
        pytest.skip("Recorded Lean Nori.")
    else:
        ref_path = os.path.join(DATA_DIR, "row75.bin")
        ref_row = np.fromfile(ref_path, dtype=np.float16)
        assert np.allclose(anchor_row, ref_row)
        
        with open(META_FILE, "r") as f:
            meta = json.load(f)
        
        assert pytest.approx(gap_scaled) == meta["nori"]["gap_scaled"]
        assert pytest.approx(mu_scaled) == meta["nori"]["mu_scaled"]

