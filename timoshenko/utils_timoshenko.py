import matplotlib.pyplot as plt
import numpy as np

def plotter_timoshenko(xx, w, phi, title_w='w(x)', title_phi='phi(x)'):
    plt.figure()
    plt.plot(xx, w, marker='o')
    plt.xlabel('x'); plt.ylabel('w(x)')
    plt.title(title_w)
    plt.grid(True)

    plt.figure()
    plt.plot(xx, phi, marker='o')
    plt.xlabel('x'); plt.ylabel('φ(x)')
    plt.title(title_phi)
    plt.grid(True)

import numpy as np

def apply_dirichlet_block(A, b, dof, value):
    """
    Impone Dirichlet sul sistema (eliminazione tramite riga/colonna).
    Nota: modifica e ritorna A,b (copia già gestita dal chiamante).
    """
    A[dof, :] = 0.0
    A[:, dof] = 0.0
    A[dof, dof] = 1.0
    b[dof] = value
    return A, b

def apply_neumann_vector(F, dof, value):
    F[dof] += value
    return F

def apply_boundary_conditions_timoshenko(K, F, bc):
    """
    K, F : matrici globali
    bc : dizionario con chiavi 'left' e/o 'right', es:
         {'left': ('dirichlet', {'w':0.0, 'phi':0.0}),
          'right': ('dirichlet', {'w':0.0, 'phi':0.0})}

    Restituisce (Kc, Fc) con le condizioni applicate.
    """
    Kc = K.copy()
    Fc = F.copy()

    # helper per ritornare i dof di un nodo
    def node_dofs(node_index):
        return 2 * node_index, 2 * node_index + 1

    nnodes = Kc.shape[0] // 2

    # supporta sia 'left' sia 'right'
    for side in ('left', 'right'):
        if side in bc and bc[side] is not None:
            typ, vals = bc[side]
            if typ.lower() == 'dirichlet':
                # scegli indice nodo
                if side == 'left':
                    node_idx = 0
                else:
                    node_idx = nnodes - 1

                dof_w, dof_phi = node_dofs(node_idx)
                # imposta w se fornito
                w_val = vals.get('w', None)
                if w_val is not None:
                    Kc, Fc = apply_dirichlet_block(Kc, Fc, dof_w, w_val)
                # imposta phi se fornito
                phi_val = vals.get('phi', None)
                if phi_val is not None:
                    Kc, Fc = apply_dirichlet_block(Kc, Fc, dof_phi, phi_val)

            elif typ.lower() == 'neumann':
                # carichi concentrati applicati ai DOF (forza sul w o momento su phi)
                # example vals {'w': value, 'phi': value}
                if side == 'left':
                    node_idx = 0
                else:
                    node_idx = nnodes - 1
                dof_w, dof_phi = node_dofs(node_idx)
                w_val = vals.get('w', None)
                phi_val = vals.get('phi', None)
                if w_val is not None:
                    Fc = apply_neumann_vector(Fc, dof_w, w_val)
                if phi_val is not None:
                    Fc = apply_neumann_vector(Fc, dof_phi, phi_val)
            else:
                raise ValueError(f"Tipo BC non riconosciuto: {typ}")

    return Kc, Fc
