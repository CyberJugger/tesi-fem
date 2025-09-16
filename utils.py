import matplotlib.pyplot as plt
from tkinter import messagebox

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