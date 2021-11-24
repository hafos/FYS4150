import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})

def exp_epsilon(T):
    """
    T in units J/kB

    Returns the epsilon expectation value for N=4 in units J
    """
    beta = 1/T # [1/J]
    return -8*np.sinh(8*beta)/(4*(np.cosh(8*beta) + 3))

def exp_magn(T):
    """
    T in units J/kB

    Returns the |m| expectation value for N=4, unitless
    """
    beta = 1/T # [1/J]
    return (2*np.exp(8*beta) + 4)/(4*(np.cosh(8*beta)+3))

def exp_Cv(T):
    """
    T in units J/kB

    Returns Cv in unit kB, for N=4
    """
    beta = 1/T # [1/J]
    E = exp_epsilon(T)*4
    E2 = 64*np.cosh(8*beta)/(np.cosh(8*beta) + 3)
    return (E2 - E**2)/(T**2*4) # [kB]

def exp_X(T):
    """
    T in units J/kB

    Returns X in unit 1/J, for N=4
    """
    beta = 1/T # [1/J]
    M = exp_magn(T)*4
    M2 = (8*np.exp(8*beta) + 8)/((np.cosh(8*beta)+3))
    return (M2 - M**2)/(T*4) # [1/J]


## Plot for 2x2 computed values:

# Load data:
file = open('simple_expectation_value_tests.dat')
file.readline()
T = []; n_c = []; eps = []; m = []; Cv = []; X = [];
for line in file:
    T_, n_c_, eps_, m_, Cv_, X_ = line.split()
    T += [float(T_)]; n_c += [int(n_c_)];
    eps += [float(eps_)]; m += [float(m_)];
    Cv += [float(Cv_)]; X += [float(X_)];

i1 = int(len(T)/3)
i2 = 2*i1
#print(n_s[i1], n_s[i1-1])

temp = np.linspace(0.1, 50.0, 1000)
fig, ax = plt.subplots(2,2, figsize=(9, 5))
ax[0,0].plot(temp, exp_epsilon(temp), label='Analytical')
ax[0,0].plot(T[i2:], eps[i2:], '+', label='Computed (10$^6$ cycles)', markersize=9)
ax[0,0].set_ylabel(r'$\left<\epsilon\right>$ [J]')
ax[0,0].set_xlabel(r'T [J/k$_B$]')

ax[0,1].plot(temp, exp_magn(temp))
ax[0,1].plot(T[i2:], m[i2:], '+', label='Computed (%s cycles)' %n_c[-1], markersize=9)
ax[0,1].set_xlabel(r'T [J/k$_B$]')
ax[0,1].set_ylabel(r'$\left<|m|\right>$ [1]')

ax[1,0].plot(temp, exp_Cv(temp))
ax[1,0].plot(T[i2:], Cv[i2:], '+', label='Computed (%s cycles)' %n_c[-1], markersize=9)
ax[1,0].set_xlabel(r'T [J/k$_B$]')
ax[1,0].set_ylabel(r'$C_V$ [k$_B$]')

ax[1,1].plot(temp, exp_X(temp))
ax[1,1].plot(T[i2:], X[i2:], '+', label='Computed (%s cycles)' %n_c[-1], markersize=9)
ax[1,1].set_xlabel(r'T [J/k$_B$]')
ax[1,1].set_ylabel(r'$\chi$ [1/J]')


ax[0,0].legend()
plt.tight_layout()
plt.savefig("analytical.pdf")
plt.show()


plt.rcParams.update({'font.size': 14})

fig, ax = plt.subplots(2,1, figsize=(7,7))

ax[0].plot(temp, exp_Cv(temp), label='Analytical')
ax[0].plot(T[i2:], Cv[i2:], 'o', label='Computed ($10^6$ cycles)')
ax[0].plot(T[i1:i2], Cv[i1:i2], 'o', label='Computed ($10^5$ cycles)')
ax[0].plot(T[:i1], Cv[:i1], 'o', label='Computed (10$^4$ cycles)')
ax[0].set_xlabel(r'T [J/k$_B$]')
ax[0].set_ylabel(r'$C_V$ [k$_B$]')
ax[0].set_xlim([1.95,2.6])
ax[0].set_ylim([0.3,0.43])

ax[1].plot(temp, exp_X(temp))
ax[1].plot(T[i2:], X[i2:], 'o', label='Computed (%.0e cycles)' %n_c[-1])
ax[1].plot(T[i1:i2], X[i1:i2], 'o', label='Computed (%.0e cycles)' %n_c[i2-1])
ax[1].plot(T[:i1], X[:i1], 'o', label='Computed (%.0e cycles)' %n_c[i1-1])
ax[1].set_xlabel(r'T [J/k$_B$]')
ax[1].set_ylabel(r'$\chi$ [1/J]')
ax[1].set_xlim([1.95,2.6])
ax[1].set_ylim([0.08,0.14])


ax[0].legend()
plt.tight_layout()
plt.savefig("analytical_zoom.pdf")
plt.show()


## Plot burn-in tests :

plt.rcParams.update({'font.size': 10})
burnin_ordered = open('burnin_test_ordered.dat')
burnin_ordered.readline()
T = []; nc_o = []; eps_o = []; m_o = [];
for line in burnin_ordered:
    T_, nc_o_, eps_o_, m_o_ = line.split()
    T += [float(T_)]
    nc_o += [float(nc_o_)]
    eps_o += [float(eps_o_)]
    m_o += [float(m_o_)]

burnin_random = open('burnin_test_random.dat')
burnin_random.readline()
nc_r = []; eps_r = []; m_r = [];
for line in burnin_random:
    T_r_, nc_r_, eps_r_, m_r_ = line.split()
    nc_r += [float(nc_r_)]
    eps_r += [float(eps_r_)]
    m_r += [float(m_r_)]

i1 = int(len(nc_o)/2)
i2 = int(len(nc_r)/2)

fig, ax = plt.subplots(2, 2, figsize=(9, 5))
ax[0, 0].plot(nc_o[1:i1], eps_o[1:i1], label = "T = %.2f J/kB" %T[0])
ax[0, 0].plot(nc_o[i1+1:], eps_o[i1+1:], label = "T = %.2f J/kB" %T[-1])
ax[0, 0].plot(nc_o[0], eps_o[0], 'ko', markersize=3, label='Initial')
ax[0, 0].legend()
ax[0, 0].set_xlabel("Number of cycles")
ax[0, 0].set_ylabel(r"$\left<\epsilon\right>$ [J]")
ax[0, 0].set_title("Ordered")

ax[1, 0].plot(nc_o[1:i1], m_o[1:i1], label = "T = %.2f J/kB" %T[0])
ax[1, 0].plot(nc_o[i1+1:], m_o[i1+1:], label = "T = %.2f J/kB" %T[-1])
ax[1, 0].plot(nc_o[0], m_o[0], 'ko', markersize=3, label = "Initial")
ax[1, 0].set_xlabel("Number of cycles")
ax[1, 0].set_ylabel(r"$\left<|m|\right>$ [1]")

ax[0, 1].plot(nc_r[1:i2], eps_r[1:i2], label = "T = %.2f J/kB" %T[0])
ax[0, 1].plot(nc_r[i2+1:], eps_r[i2+1:], label = "T = %.2f J/kB" %T[-1])
ax[0, 1].plot(nc_r[0], eps_r[0], 'ko', markersize=3, label = "Initial")
ax[0, 1].set_xlabel("Number of cycles")
ax[0, 1].set_ylabel(r"$\left<\epsilon\right>$ [J]")
ax[0, 1].set_title("Random")

ax[1, 1].plot(nc_r[1:i2], m_r[1:i2], label = "T = %.2f J/kB" %T[0])
ax[1, 1].plot(nc_r[i2+1:], m_r[i2+1:], label = "T = %.2f J/kB" %T[-1])
ax[1, 1].plot(nc_r[0], m_r[0], 'ko', markersize=3, label = "Initial")
ax[1, 1].set_xlabel("Number of cycles")
ax[1, 1].set_ylabel(r"$\left<|m|\right>$ [1]")
plt.tight_layout()
plt.savefig("burnin.pdf")
plt.show()


## Plot histogram :


P_E = open('probability_distribution.dat')
P_E.readline()
T = []; n_c = []; E = [];
for line in P_E:
    T_, n_c_, E_, = line.split()
    T += [float(T_)]
    n_c += [float(n_c_)]
    E += [float(E_)]

i1 = int(len(T)/2)

plt.rcParams.update({'font.size': 14})
#nbins1 = len(np.unique(E[0:i1], return_counts = True)[1])
#nbins1 = 11
#nbins2 = len(np.uniqe(E[i1:], return_counts = True)[1])
#nbins2 = 180
nbins1 = np.unique(E[0:i1], return_counts = True)
nbins2 = np.unique(E[i1:], return_counts = True)

fig, ax = plt.subplots(2, 1, figsize = (7, 7))
#ax[0].hist(E[0:i1], bins = nbins1, density = True)
ax[0].bar(x = nbins1[0], height = nbins1[1]/np.sum(nbins1[1]), width = 0.01)
ax[0].set_title("T = %.2f J/kB" %T[0], fontsize = 14)
ax[0].set_xlabel("$\epsilon$ [J]", fontsize = 14)
ax[0].set_ylabel("Probability distribution [1]", fontsize = 14)
#ax[1].hist(E[i1:], bins = nbins2, density = True)
ax[1].bar(x = nbins2[0], height = nbins2[1]/np.sum(nbins1[1]), width = 0.01)
ax[1].set_title("T = %.2f J/kB" %T[-1], fontsize = 14)
ax[1].set_xlabel("$\epsilon$ [J]", fontsize = 14)
ax[1].set_ylabel("Probability distribution [1]", fontsize = 14)
plt.tight_layout()
plt.savefig("probability_distribution.pdf")
plt.show()
