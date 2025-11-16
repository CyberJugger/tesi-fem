import numpy as np

# ------------------------------------------------------------------------------
# FUNZIONI DI FORMA QUADRATICHE LAGRANGE SU xi ∈ [-1, 1]
# ------------------------------------------------------------------------------
def shape_functions_quad(xi):
    N1 = 0.5 * xi * (xi - 1)
    N2 = 1.0 - xi**2
    N3 = 0.5 * xi * (xi + 1)
    return np.array([N1, N2, N3])

def shape_function_derivatives_quad(xi):
    dN1 = xi - 0.5
    dN2 = -2.0 * xi
    dN3 = xi + 0.5
    return np.array([dN1, dN2, dN3])

# ------------------------------------------------------------------------------
# ELEMENT STIFFNESS (6×6) FOR 3-NODE TIMOSHENKO ELEMENT
# ------------------------------------------------------------------------------
def timoshenko_element_quad(coords, E, I, G, A, kappa, q=0.0, m=0.0):
    """
    coords : np.array([x1, x2, x3])
    Restituisce:
        Ke (6x6), Fe (6)
    DOF per nodo: w_i, phi_i
    """

    # Gaussian points (3-point, exact for polynomials up to order 5)
    gp = np.array([ -np.sqrt(3/5), 0.0, np.sqrt(3/5) ])
    gw = np.array([  5/9,          8/9, 5/9          ])

    Ke = np.zeros((6, 6))
    Fe = np.zeros(6)

    # --------------------------------------
    # LOOP DI INTEGRAZIONE
    # --------------------------------------
    for xi, wgt in zip(gp, gw):

        # Shape functions and derivatives
        N = shape_functions_quad(xi)         # [N1,N2,N3]
        dN_dxi = shape_function_derivatives_quad(xi)

        # Mappatura isoparametrica
        x_xi = np.dot(dN_dxi, coords)        # dx/dxi
        J = x_xi                              # jacobiano 1D
        dN_dx = dN_dxi / J                    # chain rule

        # ---------------------------------------------------------
        # MATRICE B PER TRAVE DI TIMOSHENKO
        # ordine DOF: [w1,phi1, w2,phi2, w3,phi3]
        # ---------------------------------------------------------
        B = np.zeros((2, 6))

        # taglio: gamma = dw/dx - phi
        # curvatura: kappa = dphi/dx
        for i in range(3):
            B[0, 2*i]   = dN_dx[i]      # dw/dx
            B[0, 2*i+1] = -N[i]         # -phi
            B[1, 2*i+1] = dN_dx[i]      # dphi/dx

        # matrice costitutiva
        C = np.array([
            [kappa * G * A, 0],
            [0,             E * I]
        ])

        # aggiungi alla rigidezza
        Ke += (B.T @ C @ B) * J * wgt

        #carichi distribuiti
        N_full = np.zeros(6)
        for i in range(3):
            N_full[2*i] = N[i]  # solo sui w

        Fe += N_full * q * J * wgt

        M_full = np.zeros(6)
        for i in range(3):
            M_full[2*i + 1] = N[i]  # contributi su phi

        Fe += M_full * m * J * wgt

        
    return Ke, Fe


# ------------------------------------------------------------------------------
# ASSEMBLER GLOBALE PER ELEMENTI TIMOSHENKO 3-NODI QUADRATICI
# ------------------------------------------------------------------------------
def assemble_timoshenko_global(coords, connectivity, E, I, G, A, kappa, q=0, m=0):
    """
    coords        : array nodi globali (dimensione Nnodes)
    connectivity  : array (Ne x 3) con gli indici dei nodi dell’elemento
    E,I,G,A       : valori (o funzioni) delle proprietà
    kappa         : coefficiente di taglio
    q             : carico distribuito (costante)

    Restituisce:
        K (2N x 2N), F (2N)
    """

    Nnodes = len(coords)
    Ndof = 2 * Nnodes             # w, phi per nodo
    Ne = len(connectivity)

    K = np.zeros((Ndof, Ndof))
    F = np.zeros(Ndof)

    # funzione interna per gestire costanti o funzioni
    def ensure_fun(v):
        if callable(v): return v
        return lambda x: v

    E_fun = ensure_fun(E)
    I_fun = ensure_fun(I)
    G_fun = ensure_fun(G)
    A_fun = ensure_fun(A)
    q_fun = ensure_fun(q)
    m_fun = ensure_fun(m)
    # --------------------------------------
    # LOOP SUGLI ELEMENTI
    # --------------------------------------
    for e in range(Ne):

        # nodi dell’elemento
        nodes = connectivity[e]
        xe = coords[nodes]

        xm = np.mean(xe)

        Ee = E_fun(xm)
        Ie = I_fun(xm)
        Ge = G_fun(xm)
        Ae = A_fun(xm)
        qe = q_fun(xm)
        me = m_fun(xm)
        # calcolo Ke, Fe con il TUO element solver
        Ke, Fe = timoshenko_element_quad(
            coords = xe,
            E = Ee,
            I = Ie,
            G = Ge,
            A = Ae,
            kappa = kappa,
            q = qe,
            m = me
        )

        # -----------------------------
        # ASSEMBLAGGIO
        # -----------------------------
        dof = []
        for n in nodes:
            dof += [2*n, 2*n + 1]   # w_n, phi_n

        # assembly
        for a in range(6):
            Aglob = dof[a]
            F[Aglob] += Fe[a]
            for b in range(6):
                Bglob = dof[b]
                K[Aglob, Bglob] += Ke[a, b]

    return K, F
