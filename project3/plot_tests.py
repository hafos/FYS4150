import matplotlib.pyplot as plt
import numpy as np


def analytical_solution(t, x0, z0, v0, B0, V0, d, q, m):
    """
    Returns the specific analytical solution
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

def relative_error(ana_sol, num_sol):
    """
    Returns the relative error of the numerical solution
    """
    relerr = np.linalg.norm( np.abs( (ana_sol - num_sol)/ana_sol ), axis=1 )
    return relerr


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

def read_file_double(filename):
    """
    Read the test data for 2particles_RK4 from file
    """
    file = open(filename, 'r')
    B0, V0, d, q, m = file.readline().split() # Get parameters of PT and particle
    file.readline()
    t = []; x1 = []; y1 = []; z1 = []; x2 = []; y2 = []; z2 = []
    vx1 = []; vy1 = []; vz1 = []; vx2 = []; vy2 = []; vz2 = []

    for line in file:
        t_n, x1_n, y1_n, z1_n, vx1_n, vy1_n, vz1_n, x2_n, y2_n, z2_n, vx2_n, vy2_n, vz2_n = line.split()
        t += [float(t_n)]
        x1 += [float(x1_n)]; x2 += [float(x2_n)]
        y1 += [float(y1_n)]; y2 += [float(y2_n)]
        z1 += [float(z1_n)]; z2 += [float(z2_n)]
        vx1 += [float(vx1_n)]; vx2 += [float(vx2_n)]
        vy1 += [float(vy1_n)]; vy2 += [float(vy2_n)]
        vz1 += [float(vz1_n)]; vz2 += [float(vz2_n)]

    file.close()

    t = np.array(t)
    x1 = np.array(x1); y1 = np.array(y1); z1 = np.array(z1); vx1 = np.array(vx1); vy1 = np.array(vy1); vz1 = np.array(vz1)
    x2 = np.array(x2); y2 = np.array(y2); z2 = np.array(z2); vx2 = np.array(vx2); vy2 = np.array(vy2); vz2 = np.array(vz2)
    return float(B0), float(V0), float(d), float(q), float(m), t, np.array([x1,y1,z1]), np.array([vx1, vy1, vz1]), np.array([x2, y2, z2]), np.array([vx2, vy2, vz2])


"""# Open all files: They all have same B0, V0 etc.. Only change is dt
B0, V0, d, q, m, v0, t1, RKx1, RKy1, RKz1, FEx1, FEy1, FEz1 = read_file('test_results/test_RK4-FE_dt1.0e-04.dat')
B0, V0, d, q, m, v0, t2, RKx2, RKy2, RKz2, FEx2, FEy2, FEz2 = read_file('test_results/test_RK4-FE_dt5.0e-04.dat')
B0, V0, d, q, m, v0, t3, RKx3, RKy3, RKz3, FEx3, FEy3, FEz3 = read_file('test_results/test_RK4-FE_dt1.0e-03.dat')
B0, V0, d, q, m, v0, t4, RKx4, RKy4, RKz4, FEx4, FEy4, FEz4 = read_file('test_results/test_RK4-FE_dt5.0e-03.dat')
B0, V0, d, q, m, v0, t5, RKx5, RKy5, RKz5, FEx5, FEy5, FEz5 = read_file('test_results/test_RK4-FE_dt1.0e-02.dat')

# Get analytical solution:
t = np.linspace(0,100,10000)
x_ana, y_ana, z_ana = analytical_solution(t, RKx1[0], RKz1[0], v0, B0, V0, d, q, m)
wz_expected = np.sqrt(2*q*V0/(m*d**2))

## Plot the relevant stuff for the single particle tests in 9, this is in progess:

# Check if z-frequency is as expected:
fz_expected = wz_expected / (2*np.pi) # (transform from ang. freq to freq)
print('Checking frequency in z-direction:')
print('----------------------------------')
print('Expected z-frequency: %.2e' %fz_expected)
print('Approximate frequency (counting from image): %.2e \n' %(7/10))
fig, ax = plt.subplots(2,1, figsize=(8,6))
ax[0].plot(t1 , RKz1, label='RK4 method')
ax[0].plot(t, z_ana, '--', label='Analytical Solution')
ax[1].plot(t1[:len(t1)//10] , RKz1[:len(t1)//10], label='RK4 method')
ax[1].plot(t[:len(t)//10], z_ana[:len(t)//10], '--', label='Analytical Solution')
ax[0].grid(); ax[1].grid()
ax[0].set_xlabel(r't [$\mu$s]'); ax[0].set_ylabel(r'z [$\mu$m]')
ax[1].set_xlabel(r't [$\mu$s]'); ax[1].set_ylabel(r'z [$\mu$m]')
plt.legend()
ax[0].set_title(r'Single particle, h = 10$^{-4}$')
plt.tight_layout()
plt.savefig('test_results/test_frequency_z.pdf')
plt.show()"""

# Plot relative errors for all runs and methods for single particle:
B0, V0, d, q, m, t, r1, v1, r2, v2 = read_file_double('test_results/2particles_RK4_interactions_on.dat')
plt.figure()
n = 100000
plt.plot(r1[0, :n], r1[1, :n], label='1: On')
plt.plot(r2[0, :n], r2[1, :n], label='2: On')
B0, V0, d, q, m, t, r1, v1, r2, v2 = read_file_double('test_results/2particles_RK4_interactions_off.dat')
plt.plot(r1[0, :n], r1[1, :n], label='1: Off')
plt.plot(r2[0, :n], r2[1, :n], label='2: Off')
plt.legend()
plt.axis('equal')
plt.show()


"""plt.figure()
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
"""
