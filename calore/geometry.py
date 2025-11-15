import numpy as np
import matplotlib.pyplot as plt

class Geometry1D:
    def __init__(self, x_left, x_right, N=None):
        self.x_left = x_left
        self.x_right = x_right
        if N is not None:
            self.set_N(N)
    
    #allows to set a manual linspace
    def set_linspace(self, xx):
        self.xx = xx
        self.x_left = self.xx[0]
        self.x_right = self.xx[-1]
        self.N = len(self.xx)
    
    #automatically computes the linspace
    def set_N(self,N):
        self.N = N
        self.xx = np.linspace(self.x_left,self.x_right, N)
        self.dx = self.xx[1]-self.xx[0]
             

