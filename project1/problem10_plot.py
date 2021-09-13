""" Read and plot the timing data in timing_run.dat of general and specialized algorithm """
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})
# Open numerical solution output file:
file_num = open('timing_run.dat', 'r')

# Initiate nested lists to store values:
gen = [[]]; spes = [[]]; n = []
i = 0
new_run = 0

# Loop through lines:
for line in file_num:
    # Check if new run:
    if 'n:' in line:
        misc, grid_points = line.split(': ')
        n += [float(grid_points)]
        new_run = 1
    else:
        # Split each line:
        gen_str, spes_str = line.split()
        # Convert to float and add to lists:
        if new_run==1:
            gen += [[float(gen_str)]]
            spes += [[float(spes_str)]]
            i +=1
        else:
            gen[i] += [float(gen_str)]
            spes[i] += [float(spes_str)]
        new_run=0
# fixing empty first index
gen.pop(0)
spes.pop(0)
# closing read file
file_num.close()
# getting averages
gen = np.array(gen)
gen_avg = np.mean(gen, axis=1)
spes = np.array(spes)
spes_avg = np.mean(spes, axis=1)
# plotting
# print(n, gen_avg, spes_avg)
#checking ratio """
ratio = gen_avg/spes_avg
print('Ratio of the two algorithms:\n', ratio)
plt.yscale('log')
plt.xscale('log')
plt.plot(n, gen_avg, 'bo', label='general algorithm')
plt.plot(n, spes_avg, 'ro', label='specialized algorithm')
plt.legend()
plt.xlabel(r'n in $n\times n$-matrix')
plt.ylabel('time [s]')
plt.suptitle('Comparing speed at which the general algorithm')
plt.title(r'vs. specialized solves the (-1, 2, -1) $n\times n$-matrix')
plt.savefig('timing.pdf')
plt.show()
