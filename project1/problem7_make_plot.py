""" Read and plot the data in problem2_output.dat and problem7_output.dat"""
import numpy as np
import matplotlib.pyplot as plt

# Open the output file:
file_sol = open('problem2_output.dat', 'r')

# Initiate lists to store values:
x = []; u = []

# Loop through lines in file:
for line in file_sol:
    # Split each line for x and u:
    xstr, ustr = line.split()
    # Convert to float and add to lists:
    x += [float(xstr)]
    u += [float(ustr)]

# Close solution output file:
file_sol.close()

# Open numerical solution output file:
file_num = open('problem7_output.dat', 'r')

# Initiate nested lists to store values:
y = [[]]; v = [[]]
i = 0
new_run = 0

# Loop through lines:
for line in file_num:
    # Check if new run:
    if line=='\n':
        new_run = 1
    else:
        # Split each line:
        ystr, vstr = line.split()
        # Convert to float and add to lists:
        if new_run==1:
            y += [[float(ystr)]]
            v += [[float(vstr)]]
            i +=1
        new_run=0
        y[i] += [float(ystr)]
        v[i] += [float(vstr)]

# Close numerical solution output file:
file_num.close()


# Plot x and u(x):
plt.figure()
plt.plot(x, u, label='u(x)')
for i in range(len(y)):
    n = len(y[i])
    plt.plot(y[i], v[i], '--',label=f'v$_i$: n=%g'%n)
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.title('Plot of the solution u(x) \nwith numerical aproximations v_i(n)')
plt.legend()
plt.savefig('problem7_output_plot.pdf')
plt.show()

"""
Terminal> python3 problem7_make_plot.py

"""
