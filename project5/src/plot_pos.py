import time
import sys
import os
from ODESolver import RungeKutta4 as RK4
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

try:
    datafile = sys.argv[1]
except Exception:
    print "Usage: %s datafile" % sys.argv[0]
    sys.exit(1)

data = loadtxt(datafile)

n_bodies = len(data[0])/3
n_timesteps = len(data)

r = empty((n_bodies,n_timesteps,3))

print data.shape, r.shape
for i in range(n_bodies):
    r[i] = data[:,i*3:i*3+3]



fig, ax = plt.subplots()
ax = fig.gca(projection='3d')
for i in range(n_bodies):
    ax.plot(r[i,:,0], r[i,:,1], r[i,:,2], label=i)

plt.legend()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
