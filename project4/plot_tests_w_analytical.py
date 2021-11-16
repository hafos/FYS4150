import numpy as np


def exp_epsilon(T, N):
    """
    T in units J/kB
    """
    beta = 1/T # [1/J]
    return -8*np.sinh(8*beta)/(N*(np.cosh(8*beta) + 3))

def exp_magn(T, N):
    """
    T in units J/kB
    """
    beta = 1/T # [1/J]
    return (2*np.exp(8*beta) + 1)/(N*(np.cosh(8*beta)+3))

def exp_Cv(T, N):
    """
    T in units J/kB
    """
    beta = 1/T # [1/J]
    return 0

def exp_X(T, N):
    """
    T in units J/kB
    """
    beta = 1/T # [1/J]
    return 0

print(exp_epsilon(1, 4))
print(exp_magn(1, 4))
