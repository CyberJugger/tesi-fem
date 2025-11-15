import numpy as np

class Geometry1DTimoshenko:

    def __init__(self, a, b, Ne):
        self.a = a
        self.b = b
        self.Ne = Ne


        self.xx = np.linspace(a, b, 2*Ne + 1)


        self.conn = np.zeros((Ne, 3), dtype=int)

        # Ogni elemento usa nodi = [2*i, 2*i+1, 2*i+2]
        for e in range(Ne):
            self.conn[e, 0] = 2*e
            self.conn[e, 1] = 2*e + 1
            self.conn[e, 2] = 2*e + 2
