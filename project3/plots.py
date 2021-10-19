import numpy as np
import matplotlib.pyplot as plt

file = open('euler_test.dat', 'r')
t = []
x = []
y = []
z = []
x2 = []
y2 = []
z2 = []


for line in file:
    t_new, x_new, y_new, z_new, x2_new, y2_new, z2_new = line.split()
    t += [float(t_new)]
    x += [float(x_new)]
    y += [float(y_new)]
    z += [float(z_new)]
    x2 += [float(x2_new)]
    y2 += [float(y2_new)]
    z2 += [float(z2_new)]

t = np.array(t)
x = np.array(x)
y = np.array(y)
z = np.array(z)
x2 = np.array(x2)
y2 = np.array(y2)

plt.figure()
plt.plot(x,y)
plt.plot(x2,y2)
plt.axis('equal')
plt.show()
