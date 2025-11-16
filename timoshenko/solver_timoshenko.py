import numpy as np
from scipy.linalg import solve

from assembler_timoshenko import assemble_timoshenko_global   # <--- usa il TUO assembler
import utils_timoshenko as ut_t

def solver_timoshenko_1D(geom, params, bc):
    """
    geom : oggetto Geometry1DTimoshenko
        - geom.xx       = array nodi
        - geom.conn     = connettività elementi (Ne x 3)
    params :
        [E_fun, I_fun, G_fun, A_fun, kappa, q_fun]
    bc : dizionario {'left':(...), 'right':(...)}
    """

    # -------------------------
    # Estrae dati
    # -------------------------
    x = geom.xx
    conn = geom.conn

    E_fun, I_fun, G_fun, A_fun, kappa, q_fun, m_fun = params

    # -------------------------
    # ASSEMBLAGGIO GLOBALE
    # -------------------------
    K, F = assemble_timoshenko_global(
        coords = x,
        connectivity = conn,
        E = E_fun,
        I = I_fun,
        G = G_fun,
        A = A_fun,
        kappa = kappa,
        q = q_fun,
        m = m_fun
    )

    # -------------------------
    # APPLICO BOUNDARY CONDITIONS
    # -------------------------
    Kc, Fc = ut_t.apply_boundary_conditions_timoshenko(K, F, bc)

    # -------------------------
    # SOLUZIONE
    # -------------------------
    U = solve(Kc, Fc)

    # U = [w1,phi1, w2,phi2, ..., wN,phiN]
    w = U[0::2]
    phi = U[1::2]

    # -------------------------
    # PLOT
    # -------------------------
    ut_t.plotter_timoshenko(x, w, phi)

    return x, w, phi
