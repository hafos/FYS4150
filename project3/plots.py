import numpy as np
import matplotlib.pyplot as plt

def analytical_solution(t, x0, z0, v0, B0, V0, d, q, m):
    """
    The specific analytical solution
    for starting point r = (x0, 0, z0), v = (0, v0, o)
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

t = np.linspace(0, 1, 5000)
B0 = 9.65*1e1;
V0 = 9.65*1e8;
d = 1e4;
q = 1
m = 1

x0 = 1
z0 = 0
v0 = 0.1

x_ana, y_ana, z_ana = analytical_solution(t, x0, z0, v0, B0, V0, d, q, m)
plt.figure()
plt.plot(x_ana, y_ana, label='Expected')
#plt.show()

def read_file(filename):
    # Plot the two particle trajectories in the euler_test file.
    file = open(filename, 'r')
    B0, V0, d, q, m = file.readline().split() # Get parameters of PT and particles
    t = []; x = []; y = []; z = []

    for line in file:
        t_new, x_new, y_new, z_new = line.split()
        t += [float(t_new)]
        x += [float(x_new)]
        y += [float(y_new)]
        z += [float(z_new)]

    t = np.array(t)
    x = np.array(x); y = np.array(y); z = np.array(z)
    return float(B0), float(V0), float(d), float(q), float(m), t, x, y, z

B0, V0, d, q, m, t, x, y, z = read_file('euler_test.dat')
plt.plot(x,y, label='FE - Result')

B0, V0, d, q, m, t, x, y, z = read_file('RK4_test.dat')
plt.plot(x,y, '--', label='RK4 - Result')

plt.axis('equal')
plt.legend()
plt.show()

#plot the trajectories of two interacting particles (forward Euler method)

file2 = open('euler_test_2particles.dat', 'r')
B0, V0, d, q, m = file2.readline().split() # Get parameters of PT and particles
t = []
x = []; y = []; z = []
x2 = []; y2 = []; z2 = []


for line in file2:
    t_new, x_new, y_new, z_new, x2_new, y2_new, z2_new = line.split()
    t += [float(t_new)]
    x += [float(x_new)]
    y += [float(y_new)]
    z += [float(z_new)]
    x2 += [float(x2_new)]; y2 += [float(y2_new)]; z2 += [float(z2_new)]

t = np.array(t)
x = np.array(x); y = np.array(y); z = np.array(z)
x2 = np.array(x2); y2 = np.array(y2); z2 = np.array(z2)

plt.figure()
plt.plot(x,y, label='p1 - Result')
plt.plot(x2, y2, label = 'p2 - Result')
plt.axis('equal')
plt.legend()
plt.show()