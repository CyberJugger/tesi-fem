import matplotlib.pyplot as plt
import numpy as np

#prende come input le ascisse e la funzione
def plotter(xx, f, titolo_grafico, marker=False):
    plt.figure()
    if marker:
        plt.plot(xx, f, marker='o', markersize=1)
    else:
        plt.plot(xx, f)
    plt.xlabel("x"); plt.ylabel("T(x)")
    plt.title("{}".format(titolo_grafico))
    plt.grid(True)

def plot_error(x, u_fem, sol_bvp, title="Errore FEM - Analitico", show_max=True):

    # Valuta soluzione analitica sugli stessi punti della mesh FEM
    u_bvp = sol_bvp.sol(x)[0]

    # Calcola errore
    error = u_fem - u_bvp

    plotter(x, error, "Errore FEM - Analitico")
    # Plot
 
    if show_max:
        print(f"Errore massimo assoluto: {np.max(np.abs(error)):.3e}")

    plt.show()
    return error