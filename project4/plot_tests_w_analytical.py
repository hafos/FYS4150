import numpy as np
import matplotlib.pyplot as plt

def exp_epsilon(T):
    """
    T in units J/kB
    """
    beta = 1/T # [1/J]
    return -8*np.sinh(8*beta)/(4*(np.cosh(8*beta) + 3))

def exp_magn(T):
    """
    T in units J/kB
    """
    beta = 1/T # [1/J]
    return (2*np.exp(8*beta) + 1)/(4*(np.cosh(8*beta)+3))

def exp_Cv(T):
    """
    T in units J/kB

    Returns Cv in unit kB
    """
    beta = 1/T # [1/J]
    E = exp_epsilon(T)*4
    E2 = 64*np.cosh(8*beta)/(np.cosh(8*beta) + 3)
    return (E2 - E**2)/(T**2*4) # [kB]

def exp_X(T):
    """
    T in units J/kB

    Returns X in unit 1/J
    """
    beta = 1/T # [1/J]
    M = exp_magn(T)*4
    M2 = (8*np.exp(8*beta) + 2)/((np.cosh(8*beta)+3))
    return (M2 - M**2)/(T*4) # [1/J]

# Load data:
file = open('simple_expectation_value_tests.dat')
file.readline()
T = []; n_c = []; n_s = []; eps = []; m = []; Cv = []; X = [];
for line in file:
    T_, n_c_, n_s_, eps_, m_, Cv_, X_ = line.split()
    T += [float(T_)]; n_c += [int(n_c_)];
    n_s += [int(n_s_)]; eps += [float(eps_)];
    m += [float(m_)]; Cv += [float(Cv_)]; X += [float(X_)];

i1 = int(len(T)/3)
#print(n_s[i1], n_s[i1-1])

temp = np.linspace(0.1, 50.0, 1000)
fig, ax = plt.subplots(2,2)
ax[0,0].plot(temp, exp_epsilon(temp), label='Analytical')
ax[0,0].plot(T[:i1], eps[:i1], '+', label='Computed')
ax[0,0].plot(T[i1:], eps[i1:], '+', label='Computed')
ax[0,0].set_ylabel(r'$\epsilon$ [J]')
ax[0,0].set_xlabel(r'T [J/k$_B$]')

ax[0,1].plot(temp, exp_magn(temp))
ax[0,1].plot(T[:i1], m[:i1], '+', label='Computed')
ax[0,1].plot(T[i1:], m[i1:], '+', label='Computed')
ax[0,1].set_xlabel(r'T [J/k$_B$]')
ax[0,1].set_ylabel(r'$|m|$ [1]')

ax[1,0].plot(temp, exp_Cv(temp))
ax[1,0].plot(T[:i1], Cv[:i1], '+', label='Computed')
ax[1,0].plot(T[i1:], Cv[i1:], '+', label='Computed')
ax[1,0].set_xlabel(r'T [J/k$_B$]')
ax[1,0].set_ylabel(r'$C_V$ [k$_B$]')

ax[1,1].plot(temp, exp_X(temp))
ax[1,1].plot(T[:i1], X[:i1], '+', label='Computed')
ax[1,1].plot(T[i1:], X[i1:], '+', label='Computed')
ax[1,1].set_xlabel(r'T [J/k$_B$]')
ax[1,1].set_ylabel(r'$\chi$ [1/J]')


ax[0,0].legend()
plt.tight_layout()

plt.show()
#print(exp_epsilon(1))
#print(exp_magn(1))
#print(exp_Cv(1))
#print(exp_X(1))
