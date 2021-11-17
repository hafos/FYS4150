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
    return (2*np.exp(8*beta) + 4)/(4*(np.cosh(8*beta)+3))

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
    M2 = (8*np.exp(8*beta) + 8)/((np.cosh(8*beta)+3))
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
i2 = 2*i1
#print(n_s[i1], n_s[i1-1])

temp = np.linspace(0.1, 50.0, 1000)
fig, ax = plt.subplots(2,2)
ax[0,0].plot(temp, exp_epsilon(temp), label='Analytical')
ax[0,0].plot(T[:i1], eps[:i1], '+', label='Computed (%s cycles)' %n_c[i1-1])
ax[0,0].plot(T[i1:i2], eps[i1:i2], '+', label='Computed (%s cycles)' %n_c[i2-1])
ax[0,0].plot(T[i2:], eps[i2:], '+', label='Computed (%s cycles)' %n_c[-1])
ax[0,0].set_ylabel(r'$\epsilon$ [J]')
ax[0,0].set_xlabel(r'T [J/k$_B$]')

ax[0,1].plot(temp, exp_magn(temp))
ax[0,1].plot(T[:i1], m[:i1], '+', label='Computed (%s cycles)' %n_c[i1-1])
ax[0,1].plot(T[i1:i2], m[i1:i2], '+', label='Computed (%s cycles)' %n_c[i2-1])
ax[0,1].plot(T[i2:], m[i2:], '+', label='Computed (%s cycles)' %n_c[-1])
ax[0,1].set_xlabel(r'T [J/k$_B$]')
ax[0,1].set_ylabel(r'$|m|$ [1]')

ax[1,0].plot(temp, exp_Cv(temp))
ax[1,0].plot(T[:i1], Cv[:i1], '+', label='Computed (%s cycles)' %n_c[i1-1])
ax[1,0].plot(T[i1:i2], Cv[i1:i2], '+', label='Computed (%s cycles)' %n_c[i2-1])
ax[1,0].plot(T[i2:], Cv[i2:], '+', label='Computed (%s cycles)' %n_c[-1])
ax[1,0].set_xlabel(r'T [J/k$_B$]')
ax[1,0].set_ylabel(r'$C_V$ [k$_B$]')

ax[1,1].plot(temp, exp_X(temp))
ax[1,1].plot(T[:i1], X[:i1], '+', label='Computed (%s cycles)' %n_c[i1-1])
ax[1,1].plot(T[i1:i2], X[i1:i2], '+', label='Computed (%s cycles)' %n_c[i2-1])
ax[1,1].plot(T[i2:], X[i2:], '+', label='Computed (%s cycles)' %n_c[-1])
ax[1,1].set_xlabel(r'T [J/k$_B$]')
ax[1,1].set_ylabel(r'$\chi$ [1/J]')


ax[0,0].legend()
plt.tight_layout()

plt.show()
#print(exp_epsilon(1))
#print(exp_magn(1))
#print(exp_Cv(1))
#print(exp_X(1))

#plot burn-in tests
burnin_ordered = open('burnin_test_ordered.dat')
burnin_ordered.readline()
T = []; nc_o = []; eps_o = []; m_o = [];
for line in burnin_ordered:
    T_, nc_o_, eps_o_, m_o_, = line.split()
    T += [float(T_)]
    nc_o += [float(nc_o_)]
    eps_o += [float(eps_o_)]
    m_o += [float(m_o_)]

burnin_random = open('burnin_test_random.dat')
burnin_random.readline()
nc_r = []; eps_r = []; m_r = [];
for line in burnin_random:
    T_r_, nc_r_, eps_r_, m_r_, = line.split()
    nc_r += [float(nc_r_)]
    eps_r += [float(eps_r_)]
    m_r += [float(m_r_)]

i1 = int(len(T)/2)

fig, ax = plt.subplots(2, 2)
ax[0, 0].plot(nc_o[:i1], eps_o[:i1], label = "T = %.2f J/kB" %T[0])
ax[0, 0].plot(nc_o[i1:], eps_o[i1:], label = "T = %.2f J/kB" %T[-1])
ax[0, 0].legend()
ax[0, 0].set_xlabel("Number of cycles")
ax[0, 0].set_ylabel(r"$\left<\epsilon\right>$ [J]")
ax[0, 0].set_title("Ordered")

ax[1, 0].plot(nc_o[:i1], m_o[:i1], label = "T = %.2f J/kB" %T[0])
ax[1, 0].plot(nc_o[i1:], m_o[i1:], label = "T = %.2f J/kB" %T[-1])
ax[1, 0].legend()
ax[1, 0].set_xlabel("Number of cycles")
ax[1, 0].set_ylabel(r"$\left<|m|\right>$ [1]")

ax[0, 1].plot(nc_r[:i1], eps_r[:i1], label = "T = %.2f J/kB" %T[0])
ax[0, 1].plot(nc_r[i1:], eps_r[i1:], label = "T = %.2f J/kB" %T[-1])
ax[0, 1].legend()
ax[0, 1].set_xlabel("Number of cycles")
ax[0, 1].set_ylabel(r"$\left<\epsilon\right>$ [J]")
ax[0, 1].set_title("Random")

ax[1, 1].plot(nc_r[:i1], m_r[:i1], label = "T = %.2f J/kB" %T[0])
ax[1, 1].plot(nc_r[i1:], m_r[i1:], label = "T = %.2f J/kB" %T[-1])
ax[1, 1].legend()
ax[1, 1].set_xlabel("Number of cycles")
ax[1, 1].set_ylabel(r"$\left<|m|\right>$ [1]")
plt.tight_layout()
plt.show()
