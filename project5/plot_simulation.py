import numpy as np
import pyarma as pa
import math

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.animation import FuncAnimation

from animate import animate_result # Animation function example

plt.rcParams.update({'font.size': 14})

prob_7a = pa.cube()
prob_7a.load("probability_prob7a.bin")
prob_7a = np.array(prob_7a)

total_prob_7a = np.zeros(prob_7a.shape[0])
for i in range(prob_7a.shape[0]):
    total_prob_7a[i] = np.sum(prob_7a[i])

t = np.linspace(0, 0.008, prob_7a.shape[0])
y = np.linspace(0, 1, prob_7a.shape[1])
x = np.linspace(0, 1, prob_7a.shape[2])
extent = [0, 1, 0, 1]

# same thing but with double slit this time

prob_7b = pa.cube()
prob_7b.load("probability_prob7b.bin")
prob_7b = np.array(prob_7b)

total_prob_7b = np.zeros(prob_7b.shape[0])
for i in range(prob_7b.shape[0]):
    total_prob_7b[i] = np.sum(prob_7b[i])


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

# looking at the imaginary part
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
plt.show()

# looking at the real part
prob_8_real = pa.cube()
prob_8_real.load("wavefunc_prob8_real.bin")
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

# looking at the total probability function
prob_8 = prob_8_real**2 + prob_8_imag**2 # Complex conjugate

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

x08 = np.argmin(np.abs(x-0.8)) # argument position where x = 0.8
prob_slice = np.trapz(prob_8[-1][:][x08])
y = np.linspace(0, 1, len(prob_8[-1][:][x08]))
plt.plot(y, prob_8[-1][:][x08]/np.trapz(prob_8[-1][:][x08]))
plt.title('Double slit detection at t = 0.002s')
plt.savefig('2_slit_detection.pdf')
plt.show()

# print(np.trapz(prob_slice))

prob_9_1_slit = pa.cube()
prob_9_1_slit.load("probability_prob9_1_slit.bin")
prob_9_1_slit = np.array(prob_9_1_slit)

# prob_slice2 = np.trapz(prob_9_1_slit[-1][:][160])
plt.plot(y, prob_9_1_slit[-1][:][x08]/np.trapz(prob_9_1_slit[-1][:][x08]))
plt.title('Single slit detection at t = 0.002s')
plt.savefig('1_slit_detection.pdf')
plt.show()

# prob_slice3 =
prob_9_3_slits = pa.cube()
prob_9_3_slits.load("probability_prob9_3_slits.bin")
prob_9_3_slits = np.array(prob_9_3_slits)

plt.plot(y, prob_9_3_slits[-1][:][x08]/np.trapz(prob_9_3_slits[-1][:][x08 ]))
plt.title('Three slit detection at t = 0.002s')
plt.savefig('3_slit_detection.pdf')
plt.show()

animate_result(x, y, np.swapaxes(prob_7a, 1, 2), t, './animation_7a.mp4')
animate_result(x, y, np.swapaxes(prob_7b, 1, 2), t, './animation_7b.mp4')
animate_result(x, y, np.swapaxes(prob_8, 1, 2), t, './animation_8.mp4')
animate_result(x, y, np.swapaxes(prob_9_1_slit, 1, 2), t, './animation_9_1_slit.mp4')
animate_result(x, y, np.swapaxes(prob_9_3_slits, 1, 2), t, './animation_9_3_slits.mp4')
