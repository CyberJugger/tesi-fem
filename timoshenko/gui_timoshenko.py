import tkinter as tk
from tkinter import ttk
import numpy as np
import geometry_timoshenko as g_t
import solver_timoshenko as sv_t
import matplotlib.pyplot as plt

def run_solver():
    try:
        x_start = float(entry_x_start.get())
        x_end = float(entry_x_end.get())
        n = int(entry_n.get())
    except Exception as e:
        print('Errore parametri dominio:', e)
        return

    E_expr = entry_E.get()
    I_expr = entry_I.get()
    G_expr = entry_G.get()
    A_expr = entry_A.get()
    kappa_expr = entry_kappa.get()
    q_expr = entry_q.get()
    m_expr = entry_m.get()

    try:
        E_fun = lambda x: eval(E_expr, {'x': x, 'np': np})
        I_fun = lambda x: eval(I_expr, {'x': x, 'np': np})
        G_fun = lambda x: eval(G_expr, {'x': x, 'np': np})
        A_fun = lambda x: eval(A_expr, {'x': x, 'np': np})
        kappa = float(eval(kappa_expr, {'np': np}))
        q_fun = lambda x: eval(q_expr, {'x': x, 'np': np})
        m_fun = lambda x: eval(m_expr, {'x': x, 'np': np})
    except Exception as e:
        print('Errore nelle funzioni di input:', e)
        return

    geom = g_t.Geometry1DTimoshenko(x_start, x_end, n)

    left_bc = ('dirichlet', {'w': float(entry_left_w.get()) if entry_left_w.get() != '' else None,
                              'phi': float(entry_left_phi.get()) if entry_left_phi.get() != '' else None})
    right_bc = ('dirichlet', {'w': float(entry_right_w.get()) if entry_right_w.get() != '' else None,
                               'phi': float(entry_right_phi.get()) if entry_right_phi.get() != '' else None})

    params = [E_fun, I_fun, G_fun, A_fun, kappa, q_fun, m_fun]

    x, w, phi = sv_t.solver_timoshenko_1D(geom, params, {'left': left_bc, 'right': right_bc})
    plt.show()
    
root = tk.Tk()
root.title('FEM Timoshenko - GUI')

frame = ttk.LabelFrame(root, text='Dominio / Mesh')
frame.pack(padx=6, pady=6, fill='x')

ttk.Label(frame, text='x start').grid(row=0, column=0)
entry_x_start = ttk.Entry(frame); entry_x_start.insert(0, '0.0'); entry_x_start.grid(row=0, column=1)

ttk.Label(frame, text='x end').grid(row=1, column=0)
entry_x_end = ttk.Entry(frame); entry_x_end.insert(0, '1.0'); entry_x_end.grid(row=1, column=1)

ttk.Label(frame, text='n elements').grid(row=2, column=0)
entry_n = ttk.Entry(frame); entry_n.insert(0, '20'); entry_n.grid(row=2, column=1)

frame_mat = ttk.LabelFrame(root, text='Materiale / Sezione')
frame_mat.pack(padx=6, pady=6, fill='x')

ttk.Label(frame_mat, text='E(x)').grid(row=0, column=0)
entry_E = ttk.Entry(frame_mat, width=30); entry_E.insert(0, '210e9'); entry_E.grid(row=0, column=1)

ttk.Label(frame_mat, text='I(x)').grid(row=1, column=0)
entry_I = ttk.Entry(frame_mat, width=30); entry_I.insert(0, '1e-6'); entry_I.grid(row=1, column=1)

ttk.Label(frame_mat, text='G(x)').grid(row=2, column=0)
entry_G = ttk.Entry(frame_mat, width=30); entry_G.insert(0, '80e9'); entry_G.grid(row=2, column=1)

ttk.Label(frame_mat, text='A(x)').grid(row=3, column=0)
entry_A = ttk.Entry(frame_mat, width=30); entry_A.insert(0, '0.01'); entry_A.grid(row=3, column=1)

ttk.Label(frame_mat, text='kappa').grid(row=4, column=0)
entry_kappa = ttk.Entry(frame_mat, width=10); entry_kappa.insert(0, '5/6'); entry_kappa.grid(row=4, column=1)

frame_load = ttk.LabelFrame(root, text='Carico')
frame_load.pack(padx=6, pady=6, fill='x')

ttk.Label(frame_load, text='q(x)').grid(row=0, column=0)
entry_q = ttk.Entry(frame_load, width=40); entry_q.insert(0, '0.0'); entry_q.grid(row=0, column=1)

ttk.Label(frame_load, text='m(x)').grid(row=1, column=0)
entry_m = ttk.Entry(frame_load, width=40); entry_m.insert(0, '0.0'); entry_m.grid(row=1, column=1)

frame_bc = ttk.LabelFrame(root, text='Condizioni al contorno (vuoto = None)')
frame_bc.pack(padx=6, pady=6, fill='x')

ttk.Label(frame_bc, text='left w').grid(row=0, column=0)
entry_left_w = ttk.Entry(frame_bc); entry_left_w.insert(0, '0.0'); entry_left_w.grid(row=0, column=1)
ttk.Label(frame_bc, text='left phi').grid(row=0, column=2)
entry_left_phi = ttk.Entry(frame_bc); entry_left_phi.insert(0, '0.0'); entry_left_phi.grid(row=0, column=3)

ttk.Label(frame_bc, text='right w').grid(row=1, column=0)
entry_right_w = ttk.Entry(frame_bc); entry_right_w.insert(0, '0.0');entry_right_w.grid(row=1, column=1)
ttk.Label(frame_bc, text='right phi').grid(row=1, column=2)
entry_right_phi = ttk.Entry(frame_bc); entry_right_phi.insert(0, '0.0');entry_right_phi.grid(row=1, column=3)

btn = ttk.Button(root, text='Run', command=run_solver)
btn.pack(pady=6)

root.mainloop()
