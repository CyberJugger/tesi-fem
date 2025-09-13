#questo file a tempo debito ospiterà il codice dell'interfaccia grafica

import solver as sv
import geometry as g
import matplotlib.pyplot as plt
import utils as ut

geom = g.Geometry1D(2,8,100)
k= lambda x: 5-0.6*x #conduttività
a=lambda _: 0.1 #funzione area
f= lambda x: 0.03*(x-6)**4 #funzione sorgente
#kop=[10**6, 0] #parametri robin
#gop = [-1, 0]
Tp_b = 1
params=[k,a,f,None,None, a(8)*Tp_b*k(8)]

#la funzione usa le c.c di dirichlet e neumann in modo tale
#da avere un confronto anche su come ho impostato il problema
sol = sv.make_bvp_solver(geom,a, k, f, a=0, b=6, T_a=-1, Tp_b=1)

# Plot
xx=geom.xx
xx, u_fem = sv.solver_T_1D(geom, params)
error = ut.plot_error(xx, u_fem, sol, title="Errore tra FEM e analitico")
plt.show()