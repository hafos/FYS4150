import numpy as np
import pyarma as pa
import math

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm
from cycler import cycler

from animate import animate_result # Animation function example

prob_7a = pa.cube()
prob_7a.load("probability_prob7a.bin")
prob_7a = np.array(prob_7a)
#print(prob_7a.shape)

total_prob_7a = np.zeros(prob_7a.shape[0])
for i in range(prob_7a.shape[0]):
    total_prob_7a[i] = np.sum(prob_7a[i])

t = np.linspace(0, 0.008, prob_7a.shape[0])
x = np.linspace(0, 1, prob_7a.shape[1])
y = np.linspace(0, 1, prob_7a.shape[2])
extent = [0, 1, 0, 1]


plt.figure()
plt.scatter(t, np.abs(np.ones(len(total_prob_7a)) - total_prob_7a), marker='x', s=16) # abs deviation
plt.xlabel("Time [s]")
plt.ylabel("Absolute probability deviation (no wall)")
plt.savefig('probability_deviation_no_wall.png')
plt.show()

fig, ax = plt.subplots(1, 3, figsize = (15, 5))
im0 = ax[0].imshow(prob_7a[0].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7a[0]), cmap = 'plasma')
ax[0].set_title('t=0')
im1 = ax[1].imshow(prob_7a[int(prob_7a.shape[0]/2)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7a[0]), cmap = 'plasma')
ax[1].set_title('t=0.004')
im2 = ax[2].imshow(prob_7a[prob_7a.shape[0]-1].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7a[0]), cmap = 'plasma')
ax[2].set_title('t=0.008')
[fig.colorbar(imi, ax=axi, shrink=0.8) for imi, axi in zip([im0, im1, im2], ax)]
plt.tight_layout()
plt.savefig('7a_snapshots.pdf')
plt.show()

## same thing but with double slit this time

prob_7b = pa.cube()
prob_7b.load("probability_prob7b.bin")
prob_7b = np.array(prob_7b)

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

fig, ax = plt.subplots(1, 3, figsize = (15, 5))
im0 = ax[0].imshow(prob_7b[0].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7b[0]), cmap = 'plasma')
ax[0].set_title('t=0')
im1 = ax[1].imshow(prob_7b[int(prob_7b.shape[0]/2)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7b[0]), cmap = 'plasma')
ax[1].set_title('t=0.004')
im2 = ax[2].imshow(prob_7b[prob_7b.shape[0]-1].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7b[0]), cmap = 'plasma')
ax[2].set_title('t=0.008')
[fig.colorbar(imi, ax=axi, shrink=0.8) for imi, axi in zip([im0, im1, im2], ax)]
plt.tight_layout()
plt.savefig('7b_snapshots.pdf')
plt.show()

prob_8_imag = pa.cube()
prob_8_imag.load("wavefunc_prob8_imag.bin")
prob_8_imag = np.array(prob_8_imag)

fig, ax = plt.subplots(1, 3, figsize = (15, 5))
im0 = ax[0].imshow(prob_8_imag[0].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8_imag[0]), cmap = 'plasma')
ax[0].set_title('t=0')
im1 = ax[1].imshow(prob_8_imag[int(prob_8_imag.shape[0]/2)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8_imag[0]), cmap = 'plasma')
ax[1].set_title('t=0.001')
im2 = ax[2].imshow(prob_8_imag[prob_8_imag.shape[0]-1].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_8_imag[0]), cmap = 'plasma')
ax[2].set_title('t=0.002')
[fig.colorbar(imi, ax=axi, shrink=0.8) for imi, axi in zip([im0, im1, im2], ax)]
plt.suptitle(r'Im($u_{ij}$)')
plt.tight_layout()
plt.savefig('8_imag.pdf')
plt.show()

prob_8_real = pa.cube()
prob_8_real.load("wavefunc_prob8_real.bin")
prob_8_real = np.array(prob_8_real)

prob_8 = prob_8_real**2 + prob_8_imag**2 # p = u^*u



##### Plot all in one 3x3 figure: (Alternative 1)
ticks = np.linspace(0.1, 0.9, 5) # Specifies placement of xticks/yticks for closely spaced figures.
extent = [0, 1, 0, 1] # Extent of the plot
fig, ax = plt.subplots(3, 3, figsize = (15/1.3,13.5/1.3), sharey=True, sharex=True, constrained_layout=True)
def plot_three(data, row, ax, name):
    """
    Plots data to a row in ax.
    name is the name of the data (p, Re(u), Im(u)).
    The row has common vmin and vmax so a common colorbar can be appended to [row, 2]
    """
    Vmax = np.max(data[0]) # Common upper limit to scale for the three times
    im0 = ax[row,0].imshow(data[0].transpose(), extent = extent,
                            vmin = 0, vmax = Vmax, cmap = 'plasma')
    im1 = ax[row,1].imshow(data[int(data.shape[0]/2)].transpose(), extent = extent,
                            vmin = 0, vmax = Vmax, cmap = 'plasma')
    im2 = ax[row,2].imshow(data[data.shape[0]-1].transpose(), extent = extent,
                            vmin = 0, vmax = Vmax, cmap = 'plasma')
    ax[row,0].set_title(name+ '\n$t$ = 0', y=0.8, x=0.8, color='w', fontsize=14)
    ax[row,1].set_title(name+'\n$t$ = 0.001', y=0.8, x=0.75, color='w', fontsize=14)
    ax[row,2].set_title(name+'\n$t$ = 0.002', y=0.8, x=0.75, color='w', fontsize=14)
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



# print(np.trapz(prob_slice))
x08 = np.argmin(np.abs(x-0.8)) # argument position where x = 0.8

plt.figure()

## Single-slit:
prob_9_1_slit = pa.cube()
prob_9_1_slit.load("probability_prob9_1_slit.bin")
prob_9_1_slit = np.array(prob_9_1_slit)

slice = prob_9_1_slit[-1, x08] # Extract slice at screen position
plt.plot(y, slice/np.sum(slice), label='Single-slit') # Plot normalized probability
plt.xlabel('$y$')
plt.ylabel('$p(y|x=0.8; t=0.002)$')

## Double-slit :
slice = prob_8[-1, x08] # Extract slice at screen position
plt.plot(y, slice/np.sum(slice), label = 'Double-slit', color='tab:red') # Plot normalized probability

# prob_slice3 =
prob_9_3_slits = pa.cube()
prob_9_3_slits.load("probability_prob9_3_slits.bin")
prob_9_3_slits = np.array(prob_9_3_slits)

slice = prob_9_3_slits[-1, x08] # Extract slice at screen position
plt.plot(y, slice/np.sum(slice), label='Triple-slit', color='tab:olive') # Plot normalized probability


plt.legend()
plt.title('Detection probability along detector screen')
plt.tight_layout()
plt.show()
plt.savefig('detection_screen.pdf')




""" Make animations : """
animate_result(x, y, np.swapaxes(prob_7a, 1, 2), t, './animations/animation_7a.mp4')
animate_result(x, y, np.swapaxes(prob_7b, 1, 2), t, './animations/animation_7b.mp4')
animate_result(x, y, np.swapaxes(prob_8, 1, 2), t, './animations/animation_8.mp4', fps=15)
animate_result(x, y, np.swapaxes(prob_9_1_slit, 1, 2), t, './animations/animation_9_1_slit.mp4', fps=15)
animate_result(x, y, np.swapaxes(prob_9_3_slits, 1, 2), t, './animations/animation_9_3_slits.mp4', fps=15)
