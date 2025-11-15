import assembler as am
import utils as ut
import numpy as np
from scipy.integrate import solve_bvp


def solver_T_1D(geom, params,bc):
    """
    k: funzione conduttività
    a: funzione area
    f: funzione sorgente
    """
    x = geom.xx  # mesh

    A = am.stiffness_assembler_1D(x, params[0], params[1]) # assemble stiffness
    b = am.load_assembler_1D(x, params[2])    # assemble load
    A, b = ut.apply_boundary_conditions(A, b, bc)

    Pf = np.linalg.solve(A, b)      # solve linear system

    ut.plotter(x,Pf,"T function FEM", True)
    return x, Pf

def make_bvp_solver(geom, A_fun, k_fun, f_fun, a, b, left_bc, right_bc):

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
            kop_left, gop_left = left_bc[1]
            conds.append(k_fun(a)*A_fun(a)*ya[1] + kop_left*ya[0] - kop_left*gop_left)

        # --- destra ---
        if right_bc[0] == "dirichlet":
            conds.append(yb[0] - right_bc[1])
        elif right_bc[0] == "neumann":
            q = k_fun(b) * A_fun(b) * right_bc[1]
            conds.append(k_fun(b)*A_fun(b)*yb[1] - q)
        elif right_bc[0] == "robin":
            kop_right, gop_right = right_bc[1]
            conds.append(k_fun(b)*A_fun(b)*yb[1] - ( - kop_right*yb[0] + kop_right*gop_right))

        return np.array(conds)

    # Mesh e guess iniziale
    x_mesh = geom.xx
    y_guess = np.zeros((2, x_mesh.size))

    sol = solve_bvp(ode_fun, bc, x_mesh, y_guess)
    ut.plotter(x_mesh, sol.sol(x_mesh)[0], "T function analytic")
    return sol
