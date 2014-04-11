import numpy as np

def loadtable():
    """
    Loads the nonbonded potential table and returns distance, potential and force
    """
    r, V = np.loadtxt('nb01', unpack=True)

    F = np.zeros(len(r))

    F[1:-1] = [(y[i+1] - y[i-1])/2./0.018 for i in range(len(x[1:-1]))]
    F[0] = F[1]
    F[-1] = F[-2]

    return x, V, F

def leapfrog(x, v, F):
    """leapfrog algo"""

def main(u):

    x, V, F = loadtable()
