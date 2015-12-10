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
#r[body,timestep,koordinat]
eq = len(r[0,:,0])/2

#finner r etter equilibrium
r = r[:,eq:,:]
average_r_array = []
for body_i in xrange(n_bodies):
    # first, test for ejected bodies
    if r[body_i, 0,0]+r[body_i,0,1] == 0:
        continue
    average_r = 0
    for time_i in xrange(n_timesteps-eq):
        average_r += linalg.norm(r[body_i, time_i,:])
    average_r /= n_timesteps
    average_r_array.append(average_r) 

print len(average_r_array)

fig, ax = plt.subplots()
rho = average_r_array
n, bins, patches = ax.hist(rho,
                           bins=arange(int(min(rho)), 
                                       int(max(rho))+1,0.5),
                           rwidth=1)
plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
ax.set_xlabel('r')
ax.set_ylabel('rho')
plt.show()
