import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import solver as sv
import geometry as g
import utils as ut
import numpy as np

def run_solver():
    # Recupero input GUI
    x_start = float(entry_x_start.get())
    x_end = float(entry_x_end.get())
    n = int(entry_n.get())

    left_type = combo_left.get().lower()
    left_val = float(entry_left_val.get())
    right_type = combo_right.get().lower()
    right_val = float(entry_right_val.get())

    # Funzioni k(x), a(x), f(x) prese come stringhe
    k_expr = entry_k.get()
    a_expr = entry_a.get()
    f_expr = entry_f.get()

    try:
        k = lambda x: eval(k_expr, {"x": x, "np": np})
        a = lambda x: eval(a_expr, {"x": x, "np": np})
        f = lambda x: eval(f_expr, {"x": x, "np": np})
    except Exception as e:
        print("Errore nell'input delle funzioni:", e)
        return

    # Mesh
    geom = g.Geometry1D(x_start, x_end, n)
    
    # FEM
    bc = {"left": [left_type, left_val], "right": [right_type, right_val]}
    if bc["left"][0] == "dirichlet":
        left_bc = ("dirichlet", bc["left"][1], None)
    elif bc["left"][0] == "neumann":
        left_bc = ("neumann", bc["left"][1], None)

    if bc["right"][0] == "dirichlet":
        right_bc = ("dirichlet", bc["right"][1], None)
    elif bc["right"][0] == "neumann":
        right_bc = ("neumann", bc["right"][1], None)

    # Analitico
    sol = sv.make_bvp_solver(
        geom, a, k, f, a=x_start, b=x_end,
        left_bc=left_bc, right_bc=right_bc
    )    

    #checks in order to pass q instead of T'
    if bc["left"][0] == "neumann":
        bc["left"][1] = k(x_start)*a(x_start)*bc["left"][1]
    if bc["right"][0] == "neumann":
        print(bc["right"][1])
        bc["right"][1] = k(x_end)*a(x_end)*bc["right"][1]
    # dopo avere costruito k,a,f e geom
 
    params = [k, a, f, None, None, 0]  
    xx, u_fem = sv.solver_T_1D(geom, params, bc=bc)
    

    #non uso ut perché è un double plotter
    plt.figure()
    plt.plot(xx, u_fem, "o-", label="FEM")
    plt.plot(xx, sol.sol(xx)[0], "-", label="Analitico")
    plt.legend()
    plt.grid(True)
    plt.title("Soluzione FEM vs Analitico")
    
    ut.plot_error(xx, u_fem, sol)

    plt.show()

    
# GUI
root = tk.Tk()
root.title("FEM")

# Parametri dominio
frame_domain = ttk.LabelFrame(root, text="Dominio")
frame_domain.pack(padx=10, pady=10, fill="x")

ttk.Label(frame_domain, text="x start:").grid(row=0, column=0)
entry_x_start = ttk.Entry(frame_domain)
entry_x_start.insert(2, "2")
entry_x_start.grid(row=0, column=1)

ttk.Label(frame_domain, text="x end:").grid(row=1, column=0)
entry_x_end = ttk.Entry(frame_domain)
entry_x_end.insert(8, "8")
entry_x_end.grid(row=1, column=1)

ttk.Label(frame_domain, text="Numero elementi n:").grid(row=2, column=0)
entry_n = ttk.Entry(frame_domain)
entry_n.insert(0, "100")
entry_n.grid(row=2, column=1)

# Condizioni al contorno
frame_bc = ttk.LabelFrame(root, text="Condizioni al contorno")
frame_bc.pack(padx=10, pady=10, fill="x")

# Sinistra
ttk.Label(frame_bc, text="Sinistra:").grid(row=0, column=0)
combo_left = ttk.Combobox(frame_bc, values=["Dirichlet", "Neumann"])
combo_left.current(0)
combo_left.grid(row=0, column=1)
entry_left_val = ttk.Entry(frame_bc)
entry_left_val.insert(0, "-1.0")
entry_left_val.grid(row=0, column=2)

# Destra
ttk.Label(frame_bc, text="Destra:").grid(row=1, column=0)
combo_right = ttk.Combobox(frame_bc, values=["Dirichlet", "Neumann"])
combo_right.current(1)
combo_right.grid(row=1, column=1)
entry_right_val = ttk.Entry(frame_bc)
entry_right_val.insert(0, "1.0")
entry_right_val.grid(row=1, column=2)

# Funzioni
frame_fun = ttk.LabelFrame(root, text="Funzioni")
frame_fun.pack(padx=10, pady=10, fill="x")

ttk.Label(frame_fun, text="k(x):").grid(row=0, column=0)
entry_k = ttk.Entry(frame_fun, width=40)
entry_k.insert(0, "5 - 0.6*x")
entry_k.grid(row=0, column=1)

ttk.Label(frame_fun, text="a(x):").grid(row=1, column=0)
entry_a = ttk.Entry(frame_fun, width=40)
entry_a.insert(0, "0.1")
entry_a.grid(row=1, column=1)

ttk.Label(frame_fun, text="f(x):").grid(row=2, column=0)
entry_f = ttk.Entry(frame_fun, width=40)
entry_f.insert(0, "0.03*(x-6)**4")
entry_f.grid(row=2, column=1)

# Bottone run
btn = ttk.Button(root, text="Esegui", command=run_solver)
btn.pack(pady=10)

root.mainloop()
