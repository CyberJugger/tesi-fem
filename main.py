import geometry as g
import assembler as am
import numpy as np
import matplotlib.pyplot as plt

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

    plt.plot(x, Pf, marker='o')
    plt.xlabel("x")
    plt.ylabel("T(x)")
    plt.title("T function")
    plt.grid(True)
    plt.show()

    return x, Pf

geom = g.Geometry1D(0,6,100)
k= lambda x: 5-0.6*x
#k=0
a=1
f= lambda x: 0.03*(x-6)**4
#f=lambda x: 0
kop=[10^6, 0]
gop = [-1, 0]
params=[k,a,f,kop,gop]

solver_T_1D(geom, params)