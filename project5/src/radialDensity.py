import time
import sys
import os
from ODESolver import RungeKutta4 as RK4
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# try:
#     datafile = sys.argv[1]
# except Exception:
#     print "Usage: %s datafile" % sys.argv[0]
#     sys.exit(1)

fig, ax = plt.subplots()
for N in [100,150, 200, 300]:
    datafile = 'data/Nbody_position_Verlet_-9.210340_5.000000_%d_1.dat' %N
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

    r = average_r_array
    n, bins, patches = ax.hist(r,
                           bins=arange(int(min(r)), 
                                       int(max(r))+1,0.25),
                               rwidth=1, normed=True, 
                               label= 'N = %d' %N)
    plt.setp(patches, 'alpha', 0.75)



    r = sort(r)

    #n0 = max(n)
    #r0 = r[argmax(n)]
    N_ = len(r)
    
    n0 = 1.3*N_*N_
    r0 = 1.5*N_**(-1/3.0)
    
    ax.plot(r, n0/(1+(r/(r0))**4), label= 'N = %d' %N)
    ax.set_xlabel('r')
    ax.set_ylabel('P(r)')
plt.legend()

# rho = []
# rmaxrange = arange(0.01, 30, 0.5)
# for Rmax in rmaxrange:
#     number_inside = sum(r<Rmax)
#     Volume = 4*pi*Rmax**3/3 
#     rho.append(number_inside/Volume)


# fig, ax = plt.subplots()
# ax.plot(rmaxrange, rho)
# r0 = r[argmax(rho)]
# rho0 = max(rho)
# ax.plot(rmaxrange, rho0/((rmaxrange/r0)*(1+rmaxrange/r0)**2))
# ax.set_xlabel('r')
# ax.set_ylabel('rho')
plt.show()

# #0.025/(1+(rmaxrange/8.0)**4)
