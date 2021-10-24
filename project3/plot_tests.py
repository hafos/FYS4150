import matplotlib.pyplot as plt
import numpy as np


def analytical_solution(t, x0, z0, v0, B0, V0, d, q, m):
    """
    The specific analytical solution
    for starting point r = (x0, 0, z0), v = (0, v0, o)
    for a single point
    """
    w0 = q*B0/m
    wz = np.sqrt(2*q*V0/(m*d**2))
    wp = w0/2 + np.sqrt(w0**2 - 2*wz**2)/2
    wm = w0/2 - np.sqrt(w0**2 - 2*wz**2)/2
    Ap = (v0 + wm*x0)/(wm - wp)
    Am = -(v0 + wp*x0)/(wm - wp)
    x = Ap*np.cos(wp*t) + Am*np.cos(wm*t)
    y = -Ap*np.sin(wp*t) + -Am*np.sin(wm*t)
    z = z0*np.cos(wz*t)
    return x,y,z


def read_file(filename):
    """
    Read the test data for RK4+FE from file
    """
    file = open(filename, 'r')
    B0, V0, d, q, m, v0 = file.readline().split() # Get parameters of PT and particle
    t = []; x1 = []; y1 = []; z1 = []; x2 = []; y2 = []; z2 = []

    for line in file:
        t_n, x1_n, y1_n, z1_n, x2_n, y2_n, z2_n = line.split()
        t += [float(t_n)]
        x1 += [float(x1_n)]; x2 += [float(x2_n)]
        y1 += [float(y1_n)]; y2 += [float(y2_n)]
        z1 += [float(z1_n)]; z2 += [float(z2_n)]

    file.close()

    t = np.array(t)
    x1 = np.array(x1); y1 = np.array(y1); z1 = np.array(z1)
    x2 = np.array(x2); y2 = np.array(y2); z2 = np.array(z2)
    return float(B0), float(V0), float(d), float(q), float(m), float(v0), t, x1, y1, z1, x2, y2, z2



# Open all files:
B0, V0, d, q, m, v0, t1, RKx1, RKy1, RKz1, FEx1, FEy1, FEz1 = read_file('test_RK4-FE_dt1.0e-04.dat')
t = np.linspace(0,100,10000)
x_ana, y_ana, z_ana = analytical_solution(t, RKx1[0], RKz1[0], v0, B0, V0, d, q, m)
B0, V0, d, q, m, v0, t2, RKx2, RKy2, RKz2, FEx2, FEy2, FEz2 = read_file('test_RK4-FE_dt5.0e-04.dat')
B0, V0, d, q, m, v0, t3, RKx3, RKy3, RKz3, FEx3, FEy3, FEz3 = read_file('test_RK4-FE_dt1.0e-03.dat')
B0, V0, d, q, m, v0, t4, RKx4, RKy4, RKz4, FEx4, FEy4, FEz4 = read_file('test_RK4-FE_dt5.0e-03.dat')
B0, V0, d, q, m, v0, t5, RKx5, RKy5, RKz5, FEx5, FEy5, FEz5 = read_file('test_RK4-FE_dt1.0e-02.dat')


# Plot the relevant stuff for the single particle tests in 9, this is in progess:
plt.figure()
plt.plot(t1 , RKz1, label='RK4 - Result')
plt.plot(t1 , FEz1, label='FE - Result')
plt.plot(t, z_ana, 'k--', label='Analytical Solution')

plt.axis('equal')
plt.xlabel('t'); plt.ylabel('z')
plt.legend()
plt.show()

plt.figure()
plt.plot(RKx3 , RKy3, label='RK4 - Result')
#plt.plot(FEx3 , FEy3, label='FE - Result')
plt.plot(RKx4 , RKy4, label='RK4 - Result')
#plt.plot(FEx4 , FEy4, label='FE - Result')
plt.plot(RKx5 , RKy5, label='RK4 - Result')
#plt.plot(FEx5 , FEy5, label='FE - Result')
plt.plot(x_ana, y_ana, 'k--', label='Analytical Solution')
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.legend()
plt.show()
