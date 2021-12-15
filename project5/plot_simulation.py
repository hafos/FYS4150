import numpy as np
import pyarma as pa
import math

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.animation import FuncAnimation

from animate import animate_result # Modified from animation function example

plt.rcParams.update({'font.size': 14})
plt.rcParams["figure.figsize"] = (6.5,4.5) # Set a standard figsize that is used unless specified


""" Plot deviation of total probability from 1 : """
## Loading data from first simulation without walls :
prob_7a = pa.cube()
prob_7a.load("probability_prob7a.bin")
prob_7a = np.array(prob_7a)

## Compute total probability :
total_prob_7a = np.zeros(prob_7a.shape[0])
for i in range(prob_7a.shape[0]):
    total_prob_7a[i] = np.sum(prob_7a[i])

## Construct x, y, t
t = np.linspace(0, 0.008, prob_7a.shape[0])
y = np.linspace(0, 1, prob_7a.shape[1])
x = np.linspace(0, 1, prob_7a.shape[2])

## Loading data from second simulation with double-slit :
prob_7b = pa.cube()
prob_7b.load("probability_prob7b.bin")
prob_7b = np.array(prob_7b)

## Compute total probability :
total_prob_7b = np.zeros(prob_7b.shape[0])
for i in range(prob_7b.shape[0]):
    total_prob_7b[i] = np.sum(prob_7b[i])

## Plot results :
fig, ax = plt.subplots(2, 1, figsize=(7, 7))
ax[0].scatter(t, np.abs(np.ones(len(total_prob_7a)) - total_prob_7a), marker='x', s=16, alpha=0.6)
ax[0].set_title('No potential barrier')
ax[0].set_ylabel('Absolute probability deviation', labelpad=12)
ax[0].set_xlabel('Time')
ax[1].scatter(t, np.abs(np.ones(len(total_prob_7b)) - total_prob_7b), marker='x', s=16, alpha=0.6)
ax[1].set_title('Double slit')
ax[1].set_ylabel('Absolute probability deviation')
ax[1].set_xlabel('Time')
fig.tight_layout()
plt.savefig('probability_deviation_checks.pdf')
plt.show()






""" Plot colormaps of Re(u), Im(u) and p for double-slit simulation """
# Open imaginary and real results and combine to get the probability:
prob_8_imag = pa.cube()
prob_8_imag.load("wavefunc_prob8_imag.bin")
prob_8_imag = np.array(prob_8_imag)

prob_8_real = pa.cube()
prob_8_real.load("wavefunc_prob8_real.bin")
prob_8_real = np.array(prob_8_real)

prob_8 = prob_8_real**2 + prob_8_imag**2 # p = u^*u



##### Plot all in one 3x3 figure: (Alternative 1)
ticks = np.linspace(0.1, 0.9, 5) # Specifies placement of xticks/yticks for closely spaced figures.
extent = [0, 1, 0, 1] # Extent of the plot
fig, ax = plt.subplots(3, 3, figsize = (15,13.5), sharey=True, sharex=True, constrained_layout=True)
def plot_three(data, row, ax, name):
    """
    Plots data to a row in ax.
    name is the name of the data (p, Re(u), Im(u)).
    """
    Vmax = np.max(data[0])
    im0 = ax[row,0].imshow(data[0].transpose(), extent = extent,
                            vmin = 0, vmax = Vmax, cmap = 'plasma')
    im1 = ax[row,1].imshow(data[int(data.shape[0]/2)].transpose(), extent = extent,
                            vmin = 0, vmax = Vmax, cmap = 'plasma')
    im2 = ax[row,2].imshow(data[data.shape[0]-1].transpose(), extent = extent,
                            vmin = 0, vmax = Vmax, cmap = 'plasma')
    ax[row,0].set_title(name+ '\n$t$ = 0', y=0.8, x=0.8, color='w')
    ax[row,1].set_title(name+'\n$t$ = 0.001', y=0.8, x=0.75, color='w')
    ax[row,2].set_title(name+'\n$t$ = 0.002', y=0.8, x=0.75, color='w')
    return [im0, im1, im2] # To produce colorbars later
labeli = [r'$p = u^*u$', r'Re($u$)', r'Im($u$)']
## Plot each row :
ima = plot_three(prob_8, 0, ax, labeli[0])
imb = plot_three(prob_8_real, 1, ax, labeli[1])
imc = plot_three(prob_8_imag, 2, ax, labeli[2])
## Add colorbar for each row :
[fig.colorbar(imi, ax=axi, shrink=0.9, label=labeli) for imi, axi, labeli in zip([ima[2], imb[2], imc[2]], ax[:,2], labeli)]
## Set xlabel, ylabel and move ticks :
for i in range(3):
    ax[i,0].set_ylabel('y')
    ax[i,0].set_yticks(ticks)
    ax[2,i].set_xlabel('x')
    ax[2,i].set_xticks(ticks)
## Adjust spacing (do not use plt.tight_layout here)
fig.set_constrained_layout_pads(w_pad=0.1 / 144, h_pad=0.1 / 144, hspace=0.01, wspace=0.01)
plt.savefig('nine_panels.pdf')



##### Plot individual panels of 3x1: (Alternative 2)
def plot_three_panels(data, name):
    """
    Plots the data for three different times.
    Returns figure and axis objects.

    Parameters
    ----------
    data: np.array of u(x,y,t)
    """
    extent = [0, 1, 0, 1]
    fig, ax = plt.subplots(1, 3, figsize = (15,7), sharey=True, constrained_layout=True)
    im0 = ax[0].imshow(data[0].transpose(), extent = extent,
                            vmin = 0, vmax = np.amax(data[0]), cmap = 'plasma')
    ax[0].set_title('t=0')
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('y')
    ax[0].set_xticks(ticks)
    im1 = ax[1].imshow(data[int(data.shape[0]/2)].transpose(), extent = extent,
                            vmin = 0, vmax = np.amax(data[1]), cmap = 'plasma')
    ax[1].set_title('t=0.001')
    ax[1].set_xlabel('x')
    ax[1].set_xticks(ticks)
    im2 = ax[2].imshow(data[data.shape[0]-1].transpose(), extent = extent,
                            vmin = 0, vmax = np.amax(data[2]), cmap = 'plasma')
    ax[2].set_title('t=0.002')
    ax[2].set_xlabel('x')
    ax[2].set_xticks(ticks)
    fig.set_constrained_layout_pads(w_pad=0.1 / 144, h_pad=0.1 / 144, hspace=0.02, wspace=0)
    fig.colorbar(im2, ax=ax, orientation='horizontal', label=name)
    return fig, ax

fig, ax = plot_three_panels(prob_8_imag, r'Im($u$)')
plt.savefig('three_panels_imag.pdf')

fig, ax = plot_three_panels(prob_8_real, r'Re($u$)')
plt.savefig('three_panels_real.pdf')

fig, ax = plot_three_panels(prob_8, r'$p = u^*u$')
plt.savefig('three_panels_prob.pdf')


## Old version of plot:
'''fig, ax = plt.subplots(1, 3, figsize = (13, 5))
im0 = ax[0].imshow(prob_8_real[0].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8_real[0]), cmap = 'plasma')
ax[0].set_title('t=0')
im1 = ax[1].imshow(prob_8_real[int(prob_8_real.shape[0]/2)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8_real[0]), cmap = 'plasma')
ax[1].set_title('t=0.001')
im2 = ax[2].imshow(prob_8_real[prob_8_real.shape[0]-1].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8_real[0]), cmap = 'plasma')
ax[2].set_title('t=0.002')
[fig.colorbar(imi, ax=axi, shrink=0.8) for imi, axi in zip([im0, im1, im2], ax)]
plt.suptitle(r'Re($u$)')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(1, 3, figsize = (13, 5))
im0 = ax[0].imshow(prob_8[0].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8[0]), cmap = 'plasma')
ax[0].set_title('t=0')
im1 = ax[1].imshow(prob_8[int(prob_8.shape[0]/2)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8[0]), cmap = 'plasma')
ax[1].set_title('t=0.001')
im2 = ax[2].imshow(prob_8[prob_8.shape[0]-1].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8[0]), cmap = 'plasma')
ax[2].set_title('t=0.002')
[fig.colorbar(imi, ax=axi, shrink=0.8) for imi, axi in zip([im0, im1, im2], ax)]
plt.suptitle(r'$p = u^*u$')
plt.tight_layout()
plt.show()'''






""" Plotting detection probabilities along screen for single-, double- and triple-slit cases: """
x08 = np.argmin(np.abs(x-0.8)) # argument position where x = 0.8

## Double-slit :
slice = prob_8[-1, x08] # Extract slice at screen position
plt.figure()
plt.plot(y, slice/np.sum(slice)) # Plot normalized probability
plt.title('Double-slit detection probability along screen')
plt.xlabel('$y$')
plt.ylabel('$p(y|x=0.8; t=0.002)$')
plt.tight_layout()
plt.savefig('2_slit_detection.pdf')

## Single-slit:
prob_9_1_slit = pa.cube()
prob_9_1_slit.load("probability_prob9_1_slit.bin")
prob_9_1_slit = np.array(prob_9_1_slit)

slice = prob_9_1_slit[-1, x08] # Extract slice at screen position
plt.figure()
plt.plot(y, slice/np.sum(slice)) # Plot normalized probability
plt.xlabel('$y$')
plt.ylabel('$p(y|x=0.8; t=0.002)$')
plt.title('Single-slit detection probability along screen')
plt.tight_layout()
plt.savefig('1_slit_detection.pdf')

## Triple-slit:
prob_9_3_slits = pa.cube()
prob_9_3_slits.load("probability_prob9_3_slits.bin")
prob_9_3_slits = np.array(prob_9_3_slits)

slice = prob_9_3_slits[-1, x08] # Extract slice at screen position
plt.figure()
plt.plot(y, slice/np.sum(slice)) # Plot normalized probability
plt.xlabel('$y$')
plt.ylabel('$p(y|x=0.8; t=0.002)$')
plt.title('Triple-slit detection probability along screen')
plt.tight_layout()
plt.savefig('3_slit_detection.pdf')





""" Make animations : """
#animate_result(x, y, np.swapaxes(prob_7a, 1, 2), t, './animation_7a.mp4')
#animate_result(x, y, np.swapaxes(prob_7b, 1, 2), t, './animation_7b.mp4')
#animate_result(x, y, np.swapaxes(prob_8, 1, 2), t, './animation_8.mp4')
#animate_result(x, y, np.swapaxes(prob_9_1_slit, 1, 2), t, './animation_9_1_slit.mp4')
#animate_result(x, y, np.swapaxes(prob_9_3_slits, 1, 2), t, './animation_9_3_slits.mp4')
