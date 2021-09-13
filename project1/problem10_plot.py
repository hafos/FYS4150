""" Read and plot the timing data in timing_run.dat of general and specialized algorithm """
import numpy as np
import matplotlib.pyplot as plt

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



# Close numerical solution output file:
file_num.close()

print(gen[0], spes[0], n)
