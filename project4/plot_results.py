import numpy as np
import matplotlib.pyplot as plt

# Load data:
Ls = np.array([40, 60, 80, 100])
Ts = []; epss = []; ms = []; Cvs = []; Xs = [];
for L in Ls:
    file = open('results_L' + str(L) + '.dat')
    s, n_cycles = file.readline().split()
    n_cycles = int(n_cycles)
    file.readline()
    T = []; eps = []; m = []; Cv = []; X = [];
    for line in file:
        T_, eps_, m_, Cv_, X_ = line.split()
        T += [float(T_)];
        eps += [float(eps_)]; m += [float(m_)];
        Cv += [float(Cv_)]; X += [float(X_)];
    Ts += [T];
    epss += [eps]; ms += [m];
    Cvs += [Cv]; Xs += [X];
    file.close()

fig, ax = plt.subplots(2,2)
for i in range(len(Ls)):
    ax[0,0].plot(Ts[i], epss[i], label = f"L = {Ls[i]}")
ax[0,0].set_ylabel(r'$\left<\epsilon\right>$ [J]')
ax[0,0].set_xlabel(r'T [J/k$_B$]')
ax[0,0].legend()

for i in range(len(Ls)):
    ax[0,1].plot(Ts[i], ms[i], label = f"L = {Ls[i]}")
ax[0,1].set_xlabel(r'T [J/k$_B$]')
ax[0,1].set_ylabel(r'$\left<|m|\right>$ [1]')

for i in range(len(Ls)):
    ax[1,0].plot(Ts[i], Cvs[i], label = f"L = {Ls[i]}")
ax[1,0].set_xlabel(r'T [J/k$_B$]')
ax[1,0].set_ylabel(r'$C_V$ [k$_B$]')

for i in range(len(Ls)):
    ax[1,1].plot(Ts[i], Xs[i], label = f"L = {Ls[i]}")
ax[1,1].set_xlabel(r'T [J/k$_B$]')
ax[1,1].set_ylabel(r'$\chi$ [1/J]')


#ax[0,0].legend()
plt.tight_layout()
plt.savefig("results_short.pdf")
plt.show()
