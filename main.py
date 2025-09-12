import geometry as g
import assembler as am
import utils as ut
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import solve_bvp

"""
il problema nella sua forma computazionale viene cucito sulla
equazione differenziale del calore con sorgente.

non formalizza i passaggi di:
1) problema variazionale
2) discretizzazione

ma solo di:
3) problema lineare con soluzione

i passaggi specifici del problema del caloresaranno scritti su un 
file a parte
"""

def solver_T_1D(geom, params):
    #params=[k,a,f,kop,gop]
    """
    k: funzione conduttività
    a: funzione area
    f: funzione sorgente
    kop, gop: valori delle condizioni di robin
    """
    x = geom.xx  # mesh

    A = am.stiffness_assembler_1D(x, params[0], params[1],params[3])          # assemble stiffness
    b = am.load_assembler_1D(x, params[2], params[3], params[4])    # assemble load

    Pf = np.linalg.solve(A, b)      # solve linear system

    ut.plotter(x,Pf,a,"T function FEM", True)
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
    ut.plotter(x_mesh,sol.sol(x_mesh)[0],a,"T function analytic")
    return sol



geom = g.Geometry1D(2,8,100)
k= lambda x: 5-0.6*x
a=lambda _: 0.1
f= lambda x: 0.03*(x-6)**4
kop=[10**6, 0]
gop = [-1, 0]
params=[k,a,f,kop,gop]

#la funzione usa le c.c di dirichlet e neumann in modo tale
#da avere un confronto anche su come ho impostato il problema
sol = make_bvp_solver(geom,a, k, f, a=0, b=6, T_a=-1, Tp_b=0)

# Plot
xx=geom.xx
solver_T_1D(geom, params)

plt.show()