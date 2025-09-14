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

def apply_neumann(b,start:bool,value, a_fun, k_fun, x_node):
    
    q = a_fun(x_node)*k_fun(x_node) * value
    
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

def make_bvp_solver(geom, A_fun, k_fun, f_fun, a, b, T_a=None, Tp_b=None):
    def K(x): return A_fun(x)*k_fun(x)
    
    # Derivata numerica di K
    def Kprime(x, h=1e-6):
        return (K(x+h)-K(x-h))/(2*h)
    
    def ode_fun(x, y):
        dy0 = y[1]
        dy1 = -(Kprime(x)/K(x))*y[1] - f_fun(x)/K(x)
        return np.vstack((dy0, dy1))
    
    def bc(ya, yb):
        conds = []
        if T_a is not None:
            conds.append(ya[0] - T_a)
        if Tp_b is not None:
            conds.append(yb[1] - Tp_b)
        return np.array(conds)
    
    x_mesh = geom.xx
    y_guess = np.zeros((2, x_mesh.size))
    
    sol = solve_bvp(ode_fun, bc, x_mesh, y_guess)
    ut.plotter(x_mesh,sol.sol(x_mesh)[0],"T function analytic")
    return sol



