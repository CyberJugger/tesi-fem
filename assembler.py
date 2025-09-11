import numpy as np

#non mi serve
def mass_assembler_1D(x):
    n = len(x) - 1  # number of subintervals
    M = np.zeros((n+1, n+1))  # allocate mass matrix
    
    for i in range(n):  # loop over subintervals
        h = x[i+1] - x[i]  # interval length
        
        M[i, i]     += h / 3
        M[i, i+1]   += h / 6
        M[i+1, i]   += h / 6
        M[i+1, i+1] += h / 3
    
    return M

def load_assembler_1D(x, f, kop, gop):
    n = len(x) - 1
    b = np.zeros(n+1)
    
    for i in range(n):
        h = x[i+1] - x[i]
        b[i]     += f(x[i])   * h / 2
        b[i+1]   += f(x[i+1]) * h / 2
    
    b[0] = b[0] + kop[0] + gop[0]
    b[n] = b[n] + kop[1] + gop[1]

    return b

"""
def L2_projector_1D(geom, func):
    x = geom.xx  # mesh

    M = mass_assembler_1D(x)          # assemble mass
    b = load_assembler_1D(x, func)    # assemble load

    Pf = np.linalg.solve(M, b)      # solve linear system

    plt.plot(x, Pf, marker='o')
    plt.xlabel("x")
    plt.ylabel("Pf(x)")
    plt.title("L2 Projection in 1D")
    plt.grid(True)
    plt.show()

    return x, Pf
"""

def stiffness_assembler_1D(x, k, a, kop):
    n = len(x) -1
    A = np.zeros((n+1,n+1))
    
    if not callable(a):
        c_a = a
        a = lambda _: c_a

    if not callable(k):
        c_k = k
        k = lambda _: c_k

    for i in range(n):
        h = x[i+1] - x[i]#inteval amplitdue
        xm = (x[i+1]+x[i])/2 #average x on the interval
        am, km= a(xm), k(xm) #avg value of functions

        A[i, i] = A[i, i] +(am+km)/h
        A[i,i+1] = A[i,i+1] - (am+km)/h
        A[i+1,i] = A[i+1,i] - (am+km)/h
        A[i+1,i+1] = A[i+1,i+1] + (am+km)/h
    
    A[0,0] += kop[0]
    A[n,n] += kop[1]#n e non n+1 poiché 0-99, 99 è il 100esimo elemento

    return A