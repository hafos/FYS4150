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

    Returns Cv in unit kB
    """
    beta = 1/T # [1/J]
    E = exp_epsilon(T, N)*N
    E2 = 64*np.cosh(8*beta)/(np.cosh(8*beta) + 3)
    return (E2 - E**2)/(T**2*N) # [kB]

def exp_X(T, N):
    """
    T in units J/kB

    Returns X in unit 1/J
    """
    beta = 1/T # [1/J]
    M = exp_magn(T, N)*N
    M2 = (8*np.exp(8*beta) + 2)/((np.cosh(8*beta)+3))
    return (M2 - M**2)/(T*N) # [1/J]

print(exp_epsilon(1, 4))
print(exp_magn(1, 4))
print(exp_Cv(1, 4))
print(exp_X(1, 4))
