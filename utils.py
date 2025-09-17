import matplotlib.pyplot as plt
from tkinter import messagebox
import numpy as np

#prende come input le ascisse e la funzione
def plotter(xx, f, titolo_grafico, marker=False):
    plt.figure()
    if marker:
        plt.plot(xx, f, marker='o', markersize=3)
    else:
        plt.plot(xx, f)
    plt.xlabel("x"); plt.ylabel("T(x)")
    plt.title("{}".format(titolo_grafico))
    plt.grid(True)

def plot_error(x, u_fem, sol_bvp, title="Errore FEM - Analitico"):

    # Valuta soluzione analitica sugli stessi punti della mesh FEM
    u_bvp = sol_bvp.sol(x)[0]
    # Calcola errore
    error = u_fem - u_bvp
    plotter(x, error, title)
    # Plot
    plt.show()
    return error

def check_K_positive(x, k_fun, a_fun, tol=1e-12):
    for xi in x:
        Kv = a_fun(xi) * k_fun(xi)
        if Kv <= tol:
            messagebox.showerror(
                "Errore: K(x) non positiva",
                f"K(x) = a(x)*k(x) non è positiva in x={xi:.6g} (K={Kv:.6g}).\n"
                "La conduttività effettiva deve essere > 0 su tutto il dominio."
            )
            return False
    return True

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

def apply_robin(A, b, n, values):
    kop, gop = values
    
    A[n,n] += kop
    b[n] += kop*gop

    return A,b

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
    elif bc["left"][0].lower() == "robin":
        A,b = apply_robin(A, b, 0, values=bc["left"][1])

    # Estremo destro
    if bc["right"][0].lower() == "dirichlet":
        A, b = apply_dirichlet(A, b, node=n, value=bc["right"][1])
    elif bc["right"][0].lower() == "neumann":
        b = apply_neumann(b, False, bc["right"][1])
    elif bc["right"][0].lower() == "robin":
        A,b = apply_robin(A, b, 0, values=bc["right"][1])
    
    return A, b

def parse_value(entry):
    """Parsa un valore dall'entry permettendo espressioni tipo 10**6, np.pi ecc."""
    return float(eval(entry.get(), {"np": np}))


def parse_bc(side, combo, entry_val, entry_kop=None, entry_gop=None):
    """Crea la tupla boundary condition per Dirichlet, Neumann o Robin."""
    bc_type = combo.get().lower()
    try:
        if bc_type == "robin":
            kop = parse_value(entry_kop)
            gop = parse_value(entry_gop)
            return (bc_type, (kop, gop), None)
        else:
            val = parse_value(entry_val)
            return (bc_type, val, None)
    except Exception as e:
        print(f"Errore nei valori per {side}: {e}")
        return None

