import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import solver as sv
import geometry as g
import utils as ut
import numpy as np
#bollino salvatempo se le cose vanno male


def run_solver():
    # Recupero input GUI dominio
    try:
        x_start = ut.parse_value(entry_x_start)
        x_end = ut.parse_value(entry_x_end)
        n = int(eval(entry_n.get(), {"np": np}))
    except Exception as e:
        print("Errore nei parametri del dominio:", e)
        return

    # Boundary conditions
    left_bc = ut.parse_bc("sinistra", combo_left, entry_left_val, entry_left_kop, entry_left_gop)
    right_bc = ut.parse_bc("destra", combo_right, entry_right_val, entry_right_kop, entry_right_gop)

    if left_bc is None or right_bc is None:
        return

    # Funzioni k(x), a(x), f(x)
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

    if not ut.check_K_positive(geom.xx, k, a):
        return

    # FEM: costruisco dizionario bc per solver_T_1D
    bc = {"left": left_bc, "right": right_bc}

    # Conversione Neumann: passo flusso q invece che T'
    if bc["left"][0] == "neumann":
        bc["left"] = ("neumann", k(x_start) * a(x_start) * bc["left"][1], None)
    if bc["right"][0] == "neumann":
        bc["right"] = ("neumann", k(x_end) * a(x_end) * bc["right"][1], None)

    # Analitico
    sol = sv.make_bvp_solver(
        geom, a, k, f, a=x_start, b=x_end,
        left_bc=left_bc, right_bc=right_bc
    )

    # FEM
    params = [k, a, f]  
    xx, u_fem = sv.solver_T_1D(geom, params, bc=bc)
    
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

# --- Condizioni al contorno ---
frame_bc = ttk.LabelFrame(root, text="Condizioni al contorno")
frame_bc.pack(padx=10, pady=10, fill="x")

# ------------------ SINISTRA ------------------
ttk.Label(frame_bc, text="Sinistra:").grid(row=0, column=0)
combo_left = ttk.Combobox(frame_bc, values=["Dirichlet", "Neumann", "Robin"])
combo_left.current(0)
combo_left.grid(row=0, column=1)

# Entry valore
entry_left_val = ttk.Entry(frame_bc)
entry_left_val.insert(0, "-1.0")
entry_left_val.grid(row=0, column=2)

# Campi Robin (kop, gop)
label_left_kop = ttk.Label(frame_bc, text="kop:")
entry_left_kop = ttk.Entry(frame_bc)

label_left_gop = ttk.Label(frame_bc, text="gop:")
entry_left_gop = ttk.Entry(frame_bc)

def update_left_fields(event):
    if combo_left.get().lower() == "robin":
        entry_left_val.grid_forget()
        # Mostro i campi Robin affiancati
        label_left_kop.grid(row=0, column=2, padx=(5,1))
        entry_left_kop.grid(row=0, column=3, padx=(0,10))
        label_left_gop.grid(row=0, column=4, padx=(5,1))
        entry_left_gop.grid(row=0, column=5)
    else:
        label_left_kop.grid_forget()
        entry_left_kop.grid_forget()
        label_left_gop.grid_forget()
        entry_left_gop.grid_forget()
        entry_left_val.grid(row=0, column=2)

combo_left.bind("<<ComboboxSelected>>", update_left_fields)


# ------------------ DESTRA ------------------
ttk.Label(frame_bc, text="Destra:").grid(row=1, column=0)
combo_right = ttk.Combobox(frame_bc, values=["Dirichlet", "Neumann", "Robin"])
combo_right.current(1)
combo_right.grid(row=1, column=1)

# Entry valore
entry_right_val = ttk.Entry(frame_bc)
entry_right_val.insert(0, "1.0")
entry_right_val.grid(row=1, column=2)

# Campi Robin (kop, gop)
label_right_kop = ttk.Label(frame_bc, text="kop:")
entry_right_kop = ttk.Entry(frame_bc)

label_right_gop = ttk.Label(frame_bc, text="gop:")
entry_right_gop = ttk.Entry(frame_bc)

def update_right_fields(event):
    if combo_right.get().lower() == "robin":
        entry_right_val.grid_forget()
        # Mostro i campi Robin affiancati
        label_right_kop.grid(row=1, column=2, padx=(5,1))
        entry_right_kop.grid(row=1, column=3, padx=(0,10))
        label_right_gop.grid(row=1, column=4, padx=(5,1))
        entry_right_gop.grid(row=1, column=5)
    else:
        label_right_kop.grid_forget()
        entry_right_kop.grid_forget()
        label_right_gop.grid_forget()
        entry_right_gop.grid_forget()
        entry_right_val.grid(row=1, column=2)

combo_right.bind("<<ComboboxSelected>>", update_right_fields)


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
btn.pack(pady=3)

# Bottone chiudi grafici
btn_close = ttk.Button(root, text="Chiudi grafici", command=lambda: plt.close('all'))
btn_close.pack(pady=3)

root.mainloop()
