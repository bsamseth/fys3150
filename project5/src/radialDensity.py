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

fig_dummy, dummy = plt.subplots()
fig, ax = plt.subplots()
fig, dax = plt.subplots()
for N_index, N in enumerate(range(300,501,100)):
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


    r = r[:,eq:,:]
    rc = zeros((n_timesteps-eq, 3))
    for time_i in xrange(n_timesteps-eq):
        dummy_rc = zeros(3)
        number_ejected = 0
        for body_i in xrange(n_bodies):
            if r[body_i, time_i,0]+r[body_i,time_i,1] == 0:
                number_ejected += 1
                continue
                
            dummy_rc += r[body_i, time_i,:]
        dummy_rc /= (N-number_ejected)
        rc[time_i, :] = dummy_rc
    #plt.plot(sqrt(rc[:,0]**2 + rc[:,1]**2 + rc[:,2]**2))
    #plt.show()
    
    #finner r etter equilibrium
    average_r_array = []
    for body_i in xrange(n_bodies):
        # first, test for ejected bodies
        if r[body_i, 0,0]+r[body_i,0,1] == 0:
            continue
        average_r = 0
        for time_i in xrange(n_timesteps-eq):
            average_r += linalg.norm(r[body_i, time_i,:]-rc[time_i])
        average_r /= n_timesteps
        average_r_array.append(average_r) 

    print len(average_r_array)
    N_ = len(average_r_array)
    average_r_array = asarray(average_r_array)
    r = average_r_array#*N_**(1/3.0)

    n, bins, patches = dummy.hist(r,
                           bins=arange(int(min(r)), 
                                       int(max(r))+1,0.5),
                               rwidth=1, #weights = ones(N_)/N_**2,
                                  label= 'N = %d' %N,normed=True)
    plt.setp(patches, 'alpha', 1-N_index*0.33)

    binvolums = 4*pi/3 * bins**3
    dax.plot(bins[1:], n/binvolums[1:])
    ax.loglog(bins[1:]*(N**(1/3.0)), n/binvolums[1:], '+')


    r = sort(r)

    #n0 = max(n)
    #r0 = r[argmax(n)]

    

# plt.legend()

plot_r = linspace(0.1,25*500**(1/3.0),1001)
n0 =  1
r0 =  3.5
ax.plot(plot_r, n0/(1+(plot_r/(r0))**4),label=r'$%g/(1+(r/%g)^4)$' % (n0,r0))
plot_r = linspace(0.1,25,1001)
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex = True)
#dummy.plot(plot_r, n0/(1+(plot_r/(r0))**4))
dummy.set_xlabel(r'Avstand fra massesenter, [$R_0$]',size=23)
dummy.set_ylabel(r'Nummertetthet, $n(r)$',size=23)
dummy.legend(prop={'size':16})
ax.legend(prop={'size':16})
ax.set_xlabel(r'Avstand fra massesenter, [$R_0/N^{-\frac{1}{3}}$]',size=23)
ax.set_ylabel(r'Partikkeltetthet, $\Xi(r)/N^2$',size=23)
#box = dummy.get_position()
#dummy.set_position([box.x0, box.y0*1.25, box.width, box.height])





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
