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
print(prob_7a.shape)

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
plt.ylabel("Total probability (no slit)")
plt.show()

extent = [0, 1, 0, 1]

fig, ax = plt.subplots(2, 2)
im0 = ax[0, 0].imshow(prob_7a[0].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7a[0]), cmap = 'plasma')
ax[0, 0].set_title('Initial conditions')
im1 = ax[0, 1].imshow(prob_7a[int(prob_7a.shape[0]/3)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7a[0]), cmap = 'plasma')
ax[0, 1].set_title('Snapshot 1')
im2 = ax[1, 0].imshow(prob_7a[int(prob_7a.shape[0]*2/3)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7a[0]), cmap = 'plasma')
ax[1, 0].set_title('Snapshot 2')
im3 = ax[1, 1].imshow(prob_7a[prob_7a.shape[0]-1].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7a[0]), cmap = 'plasma')
ax[1, 1].set_title('Final outcome')
[fig.colorbar(imi, ax=axi, shrink=0.8, label = 'PDF') for imi, axi in zip([im0, im1, im2, im3], ax)]
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

extent = [0, 1, 0, 1]

fig, ax = plt.subplots(2, 2)
im0 = ax[0, 0].imshow(prob_7b[0].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7b[0]), cmap = 'plasma')
ax[0, 0].set_title('Initial conditions')
im1 = ax[0, 1].imshow(prob_7b[int(prob_7b.shape[0]/3)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7b[0]), cmap = 'plasma')
ax[0, 1].set_title('Snapshot 1')
im2 = ax[1, 0].imshow(prob_7b[int(prob_7b.shape[0]*2/3)].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7b[0]), cmap = 'plasma')
ax[1, 0].set_title('Snapshot 2')
im3 = ax[1, 1].imshow(prob_7b[prob_7b.shape[0]-1].transpose(), extent = extent,
                        vmin = 0, vmax = np.amax(prob_7b[0]), cmap = 'plasma')
ax[1, 1].set_title('Final outcome')
[fig.colorbar(imi, ax=axi, shrink=0.8, label = 'PDF') for imi, axi in zip([im0, im1, im2, im3], ax)]
plt.show()

animate_result(x, y, np.swapaxes(prob_7b, 1, 2), t)
