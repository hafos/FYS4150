""" Read and plot the data in problem2_output.txt """
import numpy as np
import matplotlib.pyplot as plt

# Open the output file:
file = open('problem2_output.dat', 'r')

# Initiate lists to store values:
x = []; u = []

# Loop through lines in file:
for line in file:
    # Split each line for x and u:
    xstr, ustr = line.split()
    # Convert to float and add to lists:
    x += [float(xstr)]
    u += [float(ustr)]

# Close the output file:
file.close()

# Plot x and u(x):
plt.figure()
plt.plot(x, u)
plt.xlabel('x')
plt.ylabel('u')
plt.grid()
plt.title('Plot of the solution u(x)')
plt.savefig('problem2_output_plot.pdf')
plt.show()

"""
Terminal> python3 problem2_plot_script.py

"""
