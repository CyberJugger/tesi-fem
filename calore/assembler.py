import numpy as np

#i due assembler sono specifici ad un caso in cui io abbia due condizioni di robin
def load_assembler_1D(x, f):
    n = len(x) - 1
    b = np.zeros(n+1)
    
    for i in range(n):
        h = x[i+1] - x[i]
        b[i]     += f(x[i])   * h / 2
        b[i+1]   += f(x[i+1]) * h / 2
    
    return b

def stiffness_assembler_1D(x, k, a):
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

        A[i, i] = A[i, i] +(am*km)/h
        A[i,i+1] = A[i,i+1] - (am*km)/h
        A[i+1,i] = A[i+1,i] - (am*km)/h
        A[i+1,i+1] = A[i+1,i+1] + (am*km)/h
    
    return A