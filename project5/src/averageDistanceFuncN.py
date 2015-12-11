import time
import sys
import os
from ODESolver import RungeKutta4 as RK4
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

average_distance = []

Nlist = range(50,261,10)
for N in Nlist:
    datafile = 'data/Nbody_position_Verlet_-9.210340_5.000000_%d_1.dat' %N
    data = loadtxt(datafile)

    n_bodies = len(data[0])/3
    n_timesteps = len(data)

    r = empty((n_bodies,n_timesteps,3))

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
        average_r /= n_timesteps-eq
        average_r_array.append(average_r) 
    average_average_r = sum(average_r_array)/len(average_r_array)
    average_distance.append(average_average_r)


fig, ax = plt.subplots()
ax.semilogx(Nlist, average_distance)
ax.set_xlabel('N')
ax.set_ylabel('[R]')
plt.show()

#0.025/(1+(rmaxrange/8.0)**4)
