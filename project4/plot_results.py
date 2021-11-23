import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

# Load data:
Ls = np.array([40, 60, 80, 100])
Ts = []; epss = []; ms = []; Cvs = []; Xs = [];
for L in Ls:
    file = open('higher_res_results/results_L' + str(L) + '_cycles1M.dat')
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
    T = np.array(Ts[i])
    L_fit = np.polyfit(T, epss[i], 1)
    ax[0,0].scatter(T, epss[i], s = 5, marker = "+", label = f"L = {Ls[i]}")
    ax[0,0].plot(T, L_fit[1] + L_fit[0]*T)
ax[0,0].set_ylabel(r'$\left<\epsilon\right>$ [J]')
ax[0,0].set_xlabel(r'T [J/k$_B$]')
ax[0,0].legend()

for i in range(len(Ls)):
    T = np.array(Ts[i])
    L_fit = np.polyfit(T, ms[i], 2)
    ax[0,1].scatter(T, ms[i], s = 5, marker = "+", label = f"L = {Ls[i]}")
    ax[0,1].plot(T, L_fit[2] + L_fit[1]*T + L_fit[0]*T**2)
ax[0,1].set_xlabel(r'T [J/k$_B$]')
ax[0,1].set_ylabel(r'$\left<|m|\right>$ [1]')

T_cCv = np.zeros(len(Ls)) #critical temperature
for i in range(len(Ls)):
    T = np.array(Ts[i])
    L_fit = np.polynomial.polynomial.Polynomial.fit(T, Cvs[i], 2)
    Cvs_max = np.max(L_fit(T))
    T_cCv[i] = T[np.argmax(L_fit(T))]

    ##or with smoothing splines:
    #L_fit = interpolate.UnivariateSpline(T, Cvs[i])
    #Cvs_max = np.max(L_fit(T))
    #T_cCv[i] = T[np.argmax(L_fit(T))]

    ax[1,0].scatter(T, Cvs[i], s = 5, marker = "+", label = f"L = {Ls[i]}")
    base, = ax[1,0].plot(T, L_fit(T))
    #base, = ax[1,0].plot(T, L_fit(T))
    ax[1,0].scatter(T_cCv[i], Cvs_max, color = base.get_color())
ax[1,0].set_xlabel(r'T [J/k$_B$]')
ax[1,0].set_ylabel(r'$C_V$ [k$_B$]')

T_cX = np.zeros(len(Ls))
for i in range(len(Ls)):
    T = np.array(Ts[i])
    L_fit = np.polynomial.polynomial.Polynomial.fit(T, Xs[i], 2)
    X_max = np.max(L_fit(T))
    T_cX[i] = T[np.argmax(L_fit(T))]

    #L_fit_X = interpolate.UnivariateSpline(T, Xs[i])
    #X_max = np.max(L_fit_X(T))
    #T_cX[i] = T[np.argmax(L_fit_X(T))]

    ax[1,1].scatter(T, Xs[i], s = 5, marker = "+", label = f"L = {Ls[i]}")
    base, = ax[1,1].plot(T, L_fit(T))
    #base, = ax[1,1].plot(T, L_fit_X(T))
    ax[1,1].scatter(T_cX[i], X_max, color = base.get_color())
ax[1,1].set_xlabel(r'T [J/k$_B$]')
ax[1,1].set_ylabel(r'$\chi$ [1/J]')


#ax[0,0].legend()
plt.tight_layout()
plt.savefig("results_short.pdf")
plt.show()

#print(T_cCv)

plt.figure()
plt.scatter(Ls, T_cCv)
Tc_fit = np.polynomial.polynomial.Polynomial.fit(Ls, T_cCv, 1)
plt.plot(Ls, Tc_fit(Ls), ls = '--', color = 'k')
plt.xlabel('L')
plt.ylabel(r'T$_{\rm c}$(L) [J/k$_B$]')
plt.title(r'Critical temperature (from C$_V$)')
plt.savefig('Tc_fit.pdf')
plt.show()


#find Tc at L = infinity
#first find a:
a = (Tc_fit(Ls)[0] - Tc_fit(Ls)[-1])/0.015
#then plug back into equation for each L and average the results
Tc_infty = np.zeros(len(Ls))
for i in range(len(Ls)):
    Tc_infty[i] = T_cCv[i] - a/Ls[i]
Tc_infty_avg = np.mean(Tc_infty)
print(Tc_infty_avg)
