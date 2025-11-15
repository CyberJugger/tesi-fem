import numpy as np
import geometry_timoshenko as g_t
import solver_timoshenko as sv_t

# Esempio: trave incastrata a sinistra, caricata uniformemente
# Parametri geometrici e materiali
L = 1.0                # lunghezza [m]
N = 20                 # numero elementi
E = 210e9              # [Pa]
I = 1e-6               # [m^4]
G = 80e9               # [Pa]
A = 0.01               # [m^2]
kappa = 5/6            # coeff. correzione taglio
q0 = -1000.0           # carico uniforme [N/m] (negativo = verso il basso)

# Funzioni costanti
E_fun = lambda x: E
I_fun = lambda x: I
G_fun = lambda x: G
A_fun = lambda x: A
q_fun = lambda x: q0

# Costruzione geometria
geom = g_t.Geometry1DTimoshenko(0, L, N)

# Condizioni al contorno:
# sinistra incastrata (w=0, phi=0), destra libera
bc = {
    'left': ('dirichlet', {'w': 0.0, 'phi': 0.0}),
    'right': ('dirichlet', {'w': None, 'phi': None})
}

# Parametri pacchettizzati
params = [E_fun, I_fun, G_fun, A_fun, kappa, q_fun]

# Esegui solver
x, w, phi = sv_t.solver_timoshenko_1D(geom, params, bc)

# Stampa risultati principali
print(f"Nodi: {len(x)}")
print(f"Massimo spostamento w_max = {w.min():.6e} m (a x = {x[np.argmin(w)]:.3f})")
print(f"Massima rotazione phi_max = {np.abs(phi).max():.6e} rad")
