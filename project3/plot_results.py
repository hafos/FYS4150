import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import numpy as np


def read_freq_files(filename):
    """
    Reads in info from frequency tests
    """
    file = open(filename, 'r')
    B0, V0, d, q, m, f = file.readline().split() # Get parameters of PT and particle
    file.readline() # Skip line
    wV = []; n = [];

    for line in file:
        wV_n, n_n = line.split()
        wV += [float(wV_n)]
        n += [int(n_n)]

    file.close()

    wV = np.array(wV)
    n = np.array(n)
    return float(B0), float(V0), float(d), float(q), float(m), float(f), wV, n

def frequencies_of_system(B0, V0, d, q, m):
    """
    Find the frequencies of the system of one particle
    """
    w0 = q*B0/m
    wz = np.sqrt(2*q*V0/(m*d**2))
    wp = w0/2 + np.sqrt(w0**2 - 2*wz**2)/2
    wm = w0/2 - np.sqrt(w0**2 - 2*wz**2)/2
    return w0, wz, wp, wm

# Broad view of resonance frequencies:
"""B0, V0, d, q, m, f1, wV1, n1 = read_freq_files('test_results/resonance_interactions_off_f=1.0e-01_broad.dat')
B0, V0, d, q, m, f2, wV2, n2 = read_freq_files('test_results/resonance_interactions_off_f=4.0e-01_broad.dat')
B0, V0, d, q, m, f3, wV3, n3 = read_freq_files('test_results/resonance_interactions_off_f=7.0e-01_broad.dat')

fig, ax = plt.subplots()
ax.plot(wV1, n1/100, label='f=%1.1f'%f1)
ax.plot(wV2, n2/100, label='f=%1.1f'%f2)
ax.plot(wV3, n3/100, label='f=%1.1f'%f3)
plt.legend()
plt.title('Fraction of particles left in system after 500$\mu$s')
plt.xlabel(r'$w_V$ [MHz]')
plt.ylabel(r'n/n$_{init}$ [1]')
plt.tight_layout()
plt.savefig('resonance_frequencies.pdf')
plt.show()"""


# Narrow view of resonance frequencies with interactions
B0, V0, d, q, m, f1, wV1, n1 = read_freq_files('test_results/resonance_interactions_off_f=1.0e-01_narrow.dat')
B0, V0, d, q, m, f2, wV2, n2 = read_freq_files('test_results/resonance_interactions_off_f=4.0e-01_narrow.dat')
B0, V0, d, q, m, f3, wV3, n3 = read_freq_files('test_results/resonance_interactions_off_f=7.0e-01_narrow.dat')

B0, V0, d, q, m, f1_on, wV1_on, n1_on = read_freq_files('test_results/resonance_interactions_on_f=1.0e-01_narrow.dat')
B0, V0, d, q, m, f2_on, wV2_on, n2_on = read_freq_files('test_results/resonance_interactions_on_f=4.0e-01_narrow.dat')
B0, V0, d, q, m, f3_on, wV3_on, n3_on = read_freq_files('test_results/resonance_interactions_on_f=7.0e-01_narrow.dat')

# Print the frequencies of the system
w0, wz, wp, wm = frequencies_of_system(B0, V0, d, q, m)
print(B0, V0, d, q, m)
print('Frequencies in system:')
print('w0:', w0)
print('wz:', wz)
print('wp:', wp)
print('wm:', wm)

fig, ax = plt.subplots()
ax.plot(wV1, n1/100, 'slateblue', linestyle='-', label='f=%1.1f'%f1)
ax.plot(wV1_on, n1_on/100, 'slateblue', linestyle='--')
ax.plot(wV2, n2/100, 'orange', linestyle='-', label='f=%1.1f'%f2)
ax.plot(wV2_on, n2_on/100, 'orange', linestyle='--')
ax.plot(wV3, n3/100, 'limegreen', linestyle='-', label='f=%1.1f'%f3)
ax.plot(wV3_on, n3_on/100, 'limegreen', linestyle='--')
plt.legend()
plt.title('Fraction of particles left in system after 500$\mu$s')
plt.xlabel(r'$w_V$ [MHz]')
plt.ylabel(r'n/n$_{init}$ [1]')
plt.tight_layout()
plt.savefig('resonance_frequencies_narrow.pdf')
plt.show()
