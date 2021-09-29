""" Read and plot the transformation against size of matrix datasets """
import numpy as np
import matplotlib.pyplot as plt

# plt.rcParams.update({'font.size': 14})

def plot(plot_a=True, plot_b=False, save=False):
    """ Choose which plot to generate and save  """
    if plot_a == True:
        # Open numerical solution output file:
        file_num = open('problem6_a_output.dat', 'r')
    if plot_b == True:
        file_num = open('problem6_b_output.dat', 'r')

    i = 0
    new_run = 0

    if plot_a == True:
        iterations = []; n = []
        # Loop through lines:
        for line in file_num:
            # Check if new run:
            if 'N :' in line:
                misc, grid_points = line.split(': ')
                n += [float(grid_points)]
                new_run = 1
            else:
                # Split each line:
                misc2, value = line.split(': ')
                # Convert to float and add to lists:
                if new_run==1:
                    iterations += [float(value)]
                    i +=1
                    new_run=0
    if plot_b == True:
        iterations = [[]]; n = []
        # Loop through lines:
        for line in file_num:
            # Check if new run:
            if 'N :' in line:
                misc, grid_points = line.split(': ')
                n += [float(grid_points)]
                new_run = 1
            else:
                # Split each line:
                misc2, value = line.split(': ')
                # Convert to float and add to lists:
                if new_run==1:
                    iterations += [[float(value)]]
                    i +=1
                else:
                    iterations[i] += [float(value)]
                new_run=0
        iterations.pop(0)
        n = np.array(n)
        iterations = np.array(iterations)
        iterations_avg = np.mean(iterations, axis=1)
    file_num.close()

    # print(np.polyfit(np.log(n), np.log(iterations), 1)[0]) didn't produce a good fit so experimented
    plt.xlabel(r'N in $N\times N$-matrix')
    plt.ylabel('# of similarity transformations')
    if plot_a == True:
        plt.plot(n, iterations, 'ro', markersize=4, label='h')
        plt.plot(n, np.array(n)**2.07, label=r'$N^{2.07}$')
        plt.legend()
        plt.suptitle('Scaling of transformations against size of')
        plt.title(r' $N\times N$-matrix')
    if plot_b == True:
        plt.plot(n, iterations_avg, 'ro', markersize=4)
        plt.plot(n, np.array(n)**2.033, label=r'$N^{2.033}$')
        plt.legend()
        plt.suptitle('Scaling of transformations against size of')
        plt.title(r'a dense randomly filled $N\times N$-matrix,')
    if save == True and plot_a == True:
        plt.savefig('problem6_a.pdf')
    if save == True and plot_b == True:
        plt.savefig('problem6_b.pdf')
    plt.show()

if __name__ == '__main__':
    # plot(plot_a=True, save=True)
    plot(plot_a=False, plot_b=True, save=True)
