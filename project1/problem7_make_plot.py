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
        else:
            y[i] += [float(ystr)]
            v[i] += [float(vstr)]
        new_run=0

# Close numerical solution output file:
file_num.close()

# Plot x and u(x):
'''plt.figure()
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
plt.show()'''


# Plot log of absolute error:
plt.figure()
u_i = [] # Store u for each run length
for i in range(len(y)): # Loop over each run i
    n = len(y[i])
    k_closest = [] # Store index k of x that fit y_j here
    for j in range(len(y[i])): # Loop over all y_j
        arr = np.abs(np.array(x)-y[i][j]) # Find differences of all x to this y_j
        k_closest +=[np.argmin(arr)] # Store index k of x that fits best to y_j
    u_i += [np.array(u)[k_closest]] # Add this run length u to list with all
    plt.plot(y[i], np.log10(np.abs((u_i[i]-v[i]))),label=f'v$_i$: n=%g'%n)
plt.xlabel(f'x$_i$')
plt.ylabel(f'log$_1$$_0$($\Delta_i$)')
plt.grid()
plt.title('Logarithm of the absolute error')
plt.legend()
plt.savefig('problem8_log_abserr.pdf')
plt.show()

# Plot log of relative error:
plt.figure()
relerr = [] # Store relative errors for later
for i in range(len(y)): # Loop over each run i
    n = len(y[i])
    relerr += [np.abs((u_i[i]-v[i])/u_i[i])]
    plt.plot(y[i], np.log10(relerr[i]),label=f'v$_i$: n=%g'%n)
plt.xlabel(f'x$_i$')
plt.ylabel(f'log$_1$$_0$($\epsilon_i$)')
plt.grid()
plt.title('Logarithm of the relative error')
plt.legend()
plt.savefig('problem8_log_relerr.pdf')
plt.show()

""" Print table of max relative errors: """
print('-----------------------------------------------------')
print('Table of maximal relative errors')
print('    n:     max(epsilon_i):')
print('-----------------------------------------------------')
for i in range(len(y)): # Loop over each run i
    n = len(y[i])
    print('%5g         %5.5e' %(n, np.nanmax(relerr[i])))
print('-----------------------------------------------------')


"""
Terminal> python3 problem7_make_plot.py

"""
