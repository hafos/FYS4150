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
#print(total_prob_7a)

plt.figure()
plt.plot(t, total_prob_7a)
plt.ylim(0.5, 1.5)
plt.xlabel("Time")
plt.ylabel("Total probability (no wall)")
plt.show()

extent = [0, 1, 0, 1]

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
plt.show()

## same thing but with double slit this time

prob_7b = pa.cube()
prob_7b.load("probability_prob7b.bin")
prob_7b = np.array(prob_7b)
#print(prob_7b.shape)

total_prob_7b = np.zeros(prob_7b.shape[0])
for i in range(prob_7b.shape[0]):
    total_prob_7b[i] = np.sum(prob_7b[i])

#print(total_prob_7b)

plt.figure()
plt.plot(t, total_prob_7b)
plt.ylim(0.5, 1.5)
plt.xlabel("Time")
plt.ylabel("Total probability (double slit)")
plt.show()

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
plt.show()

prob_8 = pa.cube()
prob_8.load("probability_prob8.bin")
prob_8 = np.array(prob_8)

fig, ax = plt.subplots(1, 3, figsize = (15, 5))
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
plt.suptitle(r'$p_{ij}^n = u_{ij}^{n*} u_{ij}^n$')
plt.tight_layout()
plt.show()

prob_8_imag = pa.cube()
prob_8_imag.load("probability_prob8_imag.bin")
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
plt.show()

prob_8_real = pa.cube()
prob_8_real.load("probability_prob8_real.bin")
prob_8_real = np.array(prob_8_real)

fig, ax = plt.subplots(1, 3, figsize = (15, 5))
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
plt.suptitle(r'Re($u_{ij}$)')
plt.tight_layout()
plt.show()

#animate_result(x, y, np.swapaxes(prob_7a, 1, 2), t, './animation_7a.mp4')
#animate_result(x, y, np.swapaxes(prob_7b, 1, 2), t, './animation_7b.mp4')
animate_result(x, y, np.swapaxes(prob_8, 1, 2), t, './animation_8.mp4')
