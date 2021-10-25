import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', size=14)


def analytical_solution(t, x0, z0, v0, B0, V0, d, q, m):
    """
    Returns the specific analytical solution
    for starting point r = (x0, 0, z0), v = (0, v0, 0)
    for a single point
    """
    omega0 = q*B0/m
    omega_z = np.sqrt(2*q*V0/(m*d**2))
    wp = omega0/2 + np.sqrt(omega0**2 - 2*omega_z**2)/2
    wm = omega0/2 - np.sqrt(omega0**2 - 2*omega_z**2)/2
    Ap = (v0 + wm*x0)/(wm - wp)
    Am = -(v0 + wp*x0)/(wm - wp)
    x = Ap*np.cos(wp*t) + Am*np.cos(wm*t)
    y = -Ap*np.sin(wp*t) + -Am*np.sin(wm*t)
    z = z0*np.cos(omega_z*t)
    return x,y,z

def relative_error_dist(ana_sol, num_sol):
    """
    Returns the relative error and maximum error for the error convergence rate of the
    numerical solution compared to the analytical solution of the particle paths
    """
    diff = ana_sol - num_sol
    dist_diff = np.linalg.norm(diff, axis=0)
    ana_sol_dist = np.linalg.norm(ana_sol, axis=0)
    relerr = dist_diff/ana_sol_dist
    delta_max_k = np.max(dist_diff)
    return relerr, delta_max_k


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

if __name__ == "__main__":
    """ Problem 9: Single Particle experiments: """
    """ Getting data from single particle simulation run with different time steps and
        using either the Runge-Kutta4 or Forward-Euler integrator. """
    # B0, V0, d, q, m, v0, t1, RKx1, RKy1, RKz1, FEx1, FEy1, FEz1 = read_file('test_results/test_RK4-FE_dt1.0e-01.dat')
    # RK1 = np.array([RKx1, RKy1, RKz1]); FE1 = np.array([FEx1, FEy1, FEz1])
    # B0, V0, d, q, m, v0, t2, RKx2, RKy2, RKz2, FEx2, FEy2, FEz2 = read_file('test_results/test_RK4-FE_dt5.0e-02.dat')
    # RK2 = np.array([RKx2, RKy2, RKz2]); FE2 = np.array([FEx2, FEy2, FEz2])
    # B0, V0, d, q, m, v0, t3, RKx3, RKy3, RKz3, FEx3, FEy3, FEz3 = read_file('test_results/test_RK4-FE_dt1.0e-02.dat')
    # RK3 = np.array([RKx3, RKy3, RKz3]); FE3 = np.array([FEx3, FEy3, FEz3])
    # B0, V0, d, q, m, v0, t4, RKx4, RKy4, RKz4, FEx4, FEy4, FEz4 = read_file('test_results/test_RK4-FE_dt5.0e-03.dat')
    # RK4 = np.array([RKx4, RKy4, RKz4]); FE4 = np.array([FEx4, FEy4, FEz4])
    # B0, V0, d, q, m, v0, t5, RKx5, RKy5, RKz5, FEx5, FEy5, FEz5 = read_file('test_results/test_RK4-FE_dt1.0e-03.dat')
    # RK5 = np.array([RKx5, RKy5, RKz5]); FE5 = np.array([FEx5, FEy5, FEz5])
    # RK = [RK1, RK2, RK3, RK4, RK5]
    # FE = [FE1, FE2, FE3, FE4, FE5]
    # h_list = [0.1, 0.05, 0.01, 0.005, 0.001]
    # t_list = [t1, t2, t3, t4, t5]
    #
    # t = np.linspace(0, 100, len(RKx1))
    # omega_z_expected = np.sqrt(2*q*V0/(m*d**2))

    """ Checking if movement in z-direction is as expected given the z-frequency omega_z: """
    # fz_expected = omega_z_expected / (2*np.pi) # (transform from ang.freq to freq)
    # # We use the simulation run with the highest accuracy; RK4 with smallest time step
    # z_ana = analytical_solution(t_list[-1], RK[-1][0][0], RK[-1][2][0], v0, B0, V0, d, q, m)[2]
    # print('Checking frequency in z-direction over 100 microseconds:')
    # print('----------------------------------')
    # print(f'Expected z-frequency: {fz_expected}')
    # tol = 6e-7
    # max = np.argwhere(np.abs(RKz5 - 10) < tol) # finds all the peaks for infering the frequency
    # print(f'Approximate frequency from simulation: {1/(t5[max[-1]] - t5[max[-2]])}')
    # print(f'Difference: {fz_expected - 1/(t5[max[-1]] - t5[max[-2]])}')
    # plt.plot(t5, z_ana, label='Analytical Solution')
    # plt.plot(t5[::100] , RKz5[::100], '--', label=f'RK4 with h = {h_list[-1]}')
    # # If one wants to check the peaks are hit
    # # for top in max:
    # #     plt.plot(t5[top] , RKz5[top], 'o')
    # # plt.plot(t5[::100] , FEz5[::100], '.', markersize=2, label='FE method')
    # plt.grid();
    # plt.xlabel(r't [$\mu s$]'); plt.ylabel(r'z [$\mu m$]')
    # plt.legend()
    # plt.title(r'Oscillation in z-direction of a single $Ca^+$ atom')
    # plt.tight_layout()
    # plt.savefig('test_results/test_frequency_z.pdf')
    # plt.show()

    """ The relative errors and error convergence rate for the different integrators RK4 vs. FE """
    # relative_errors = [[],[]]; max_error = [[],[]]
    # for i in range(len(RK)):
    #     x_ana, y_ana, z_ana = analytical_solution(t_list[i], RK[i][0][0], RK[i][2][0], v0, B0, V0, d, q, m)
    #     ana_coord_RK = np.array([x_ana, y_ana, z_ana])
    #     x_ana2, y_ana2, z_ana2 = analytical_solution(t_list[i], FE[i][0][0], FE[i][2][0], v0, B0, V0, d, q, m)
    #     ana_coord_FE = np.array([x_ana2, y_ana2, z_ana2])
    #     relerr_RK, max_err_RK = relative_error_dist(ana_coord_RK, RK[i])
    #     relerr_FE, max_err_FE = relative_error_dist(ana_coord_FE, FE[i])
    #     relative_errors[0].append(relerr_RK); max_error[0].append(max_err_RK)
    #     relative_errors[1].append(relerr_FE); max_error[1].append(max_err_FE)
    #
    # convergence_rate_RK = 0
    # convergence_rate_FE = 0
    # for k in range(1, 5):
    #     convergence_rate_RK += np.log(max_error[0][k]/max_error[0][k-1])/np.log(h_list[k]/h_list[k-1])
    #     convergence_rate_FE += np.log(max_error[1][k]/max_error[1][k-1])/np.log(h_list[k]/h_list[k-1])
    # print(convergence_rate_RK*1/4)
    # print(convergence_rate_FE*1/4)
    # # print(f'Convergence rate for the two integrators when varrying step length h from {h_list[0]} to {h_list[-1]}:\n{'RungeKutta4':15s}:{'Forward-Euler':15s}:\n{(convergence_rate_RK*1/4):15.2f}{(convergence_rate_FE*1/4):15.2f}')
    # # Plot relative errors for all runs and methods for single particle:
    # fig, ax = plt.subplots(1, 2, figsize=(14/1.4, 7/1.4))#,figsize=(14/1.4, 7/1.4))
    # for i in range(5):
    #     # print(len(relative_errors[0][i]))
    #     ax[0].plot(t_list[i], relative_errors[0][i][:], label=f'h: {h_list[i]}')
    #     ax[0].set_title('Runge-Kutta4')
    #     ax[1].plot(t_list[i], relative_errors[1][i][:], label=f'h: {h_list[i]}')
    #     ax[1].set_title('Forward-Euler')
    #     ax[0].set_yscale('log'); ax[1].set_yscale('log')
    #     # # # plt.plot(x_ana, y_ana, 'k', label='Analytical Solution')
    #     # # plt.axis('equal')
    # fig.supylabel(r'Relative Error $[\log_{10}]$')
    # fig.supxlabel('Time 'r'$[\mu s]$')
    # fig.tight_layout()
    # plt.legend()
    # plt.savefig('test_results/relative_error.pdf')
    # plt.show()

    """ Problem 9: Test with Two particles with and without Coulomb interactions """
    pos1, vel1, pos2, vel2 = read_file_double('test_results/2particles_RK4_interactions_off.dat')[6:10]
    pos3, vel3, pos4, vel4 = read_file_double('test_results/2particles_RK4_interactions_on.dat')[6:10]
    fig = plt.figure(figsize = (7, 7))
    ax = plt.axes(projection = '3d')
    ax.plot3D(pos1[0], pos1[1], pos1[2], label = 'Particle 1')
    ax.plot3D(pos2[0], pos2[1], pos2[2], label = 'Particle 2')
    ax.set_xlabel('X [μm]')
    ax.set_ylabel('Y [μm]')
    ax.set_zlabel('Z [μm]')
    ax.set_title('No interaction')
    plt.legend()
    plt.savefig('test_results/2particles_3D_off.pdf')
    plt.show()

    fig = plt.figure(figsize = (7, 7))
    ax = plt.axes(projection = '3d')
    ax.plot3D(pos3[0], pos3[1], pos3[2], label = 'Particle 1')
    ax.plot3D(pos4[0], pos4[1], pos4[2], label = 'Particle 2')
    ax.set_xlabel('X [μm]')
    ax.set_ylabel('Y [μm]')
    ax.set_zlabel('Z [μm]')
    ax.set_title('Interaction')
    plt.legend()
    plt.savefig('test_results/2particles_3D_on.pdf')
    plt.show()

    '''
    plt.figure()
    # n = 100000
    plt.plot(pos1[0], pos1[1], label='Particle 1: Off')
    plt.plot(pos2[0], pos2[1], label='Particle 2: Off')
    plt.plot(pos3[0], pos3[1], label='Particle 1: On')
    plt.plot(pos4[0], pos4[1], label='Particle 2: On')
    plt.legend()
    plt.axis('equal')
    plt.show()
    '''
