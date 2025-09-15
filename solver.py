import assembler as am
import utils as ut
import numpy as np
from scipy.integrate import solve_bvp

def apply_dirichlet(A, b, node, value):
    # azzero riga
    A[node, :] = 0
    # imposto 1 sulla diagonale
    A[node, node] = 1.0
    # imposto il valore nel vettore dei carichi
    b[node] = value
    return A, b

def apply_neumann(b,start:bool,value):    
    if start:
        b[0] -= value
    else:
        b[-1] += value
    return b

def apply_boundary_conditions(A, b, bc):
    """
    Applica automaticamente le condizioni al contorno a seconda del tipo.
    bc è un dict del tipo:
    bc = {"left": ("dirichlet", value), "right": ("neumann", value)}
    """
    n = len(b) - 1
    
    # Estremo sinistro
    if bc["left"][0].lower() == "dirichlet":
        A, b = apply_dirichlet(A, b, node=0, value=bc["left"][1])
    elif bc["left"][0].lower() == "neumann":
        b = apply_neumann(b, True, bc["left"][1])
    
    # Estremo destro
    if bc["right"][0].lower() == "dirichlet":
        A, b = apply_dirichlet(A, b, node=n, value=bc["right"][1])
    elif bc["right"][0].lower() == "neumann":
        b = apply_neumann(b, False, bc["right"][1])
    
    return A, b


def solver_T_1D(geom, params,bc):
    #params=[k,a,f,kop,gop,q]
    """
    k: funzione conduttività
    a: funzione area
    f: funzione sorgente
    kop, gop: valori delle condizioni di robin
    q: calore ricavato da neumann
    """
    x = geom.xx  # mesh

    A = am.stiffness_assembler_1D(x, params[0], params[1],params[3]) # assemble stiffness
    b = am.load_assembler_1D(x, params[2], params[3], params[4])    # assemble load
    
    #A, b = apply_dirichlet(A, b, node=0, value=-1.0)
    #b = apply_neumann(b, False, params[5])
    #messo qui come test

    A, b = apply_boundary_conditions(A, b, bc)

    Pf = np.linalg.solve(A, b)      # solve linear system

    ut.plotter(x,Pf,"T function FEM", True)
    return x, Pf

def make_bvp_solver(geom, A_fun, k_fun, f_fun, a, b,
                    left_bc=("dirichlet", None, None),
                    right_bc=("dirichlet", None, None)):

    def K(x): return A_fun(x) * k_fun(x)

    # Derivata numerica di K
    def Kprime(x, h=1e-6):
        return (K(x+h) - K(x-h)) / (2*h)

    # Equazioni del sistema
    def ode_fun(x, y):
        dy0 = y[1]
        dy1 = -(Kprime(x)/K(x)) * y[1] - f_fun(x)/K(x)
        return np.vstack((dy0, dy1))

    # Condizioni al contorno
    def bc(ya, yb):
        conds = []

        # --- sinistra ---
        if left_bc[0] == "dirichlet":
            conds.append(ya[0] - left_bc[1])
        elif left_bc[0] == "neumann":
            q = k_fun(a) * A_fun(a) * left_bc[1]
            conds.append(k_fun(a)*A_fun(a)*ya[1] - q)
        elif left_bc[0] == "robin":
            alpha, g = left_bc[1], left_bc[2]
            conds.append(alpha*ya[0] + k_fun(a)*A_fun(a)*ya[1] - g)

        # --- destra ---
        if right_bc[0] == "dirichlet":
            conds.append(yb[0] - right_bc[1])
        elif right_bc[0] == "neumann":
            q = k_fun(b) * A_fun(b) * right_bc[1]
            conds.append(k_fun(b)*A_fun(b)*yb[1] - q)
        elif right_bc[0] == "robin":
            alpha, g = right_bc[1], right_bc[2]
            conds.append(alpha*yb[0] + k_fun(b)*A_fun(b)*yb[1] - g)

        return np.array(conds)

    # Mesh e guess iniziale
    x_mesh = geom.xx
    y_guess = np.zeros((2, x_mesh.size))

    sol = solve_bvp(ode_fun, bc, x_mesh, y_guess)
    ut.plotter(x_mesh, sol.sol(x_mesh)[0], "T function analytic")
    return sol
