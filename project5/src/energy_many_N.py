import time
import sys
import os
from ODESolver import RungeKutta4 as RK4
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex = True)
fig, ax = plt.subplots()
ax.hold(True)
Nlist = range(50,301,30)
for N in Nlist:
    datafile = "data/Nbody_energy_Verlet_-9.210340_5.000000_%d_1.dat" % N
    logh, Tmax, N, check_ =  [float(datafile[:-4].split('_')[i+3]) for i in range(0,4)]
    
    if check_:
        s = 'Her fjerner vi utskytne partikler'
    else:
        s = 'Her fjerner vi \_ikke\_ utskytne partikler'

        
    data = loadtxt(datafile)

    t_array = [Tmax*i/1000.0 for i in range(len(data[:,0]))]

    ax.plot(t_array, (-data[:,0]/N)/max(data[:,0]/N), label=r"$N=%d$" % N)


ax.set_title(r'Energi som funksjon av tiden', size=23)
ax.set_xlabel(r'Tid  [$\tau_{crunch}$] \\ \\' \
           + r'Parametre: $h = %g, T_{max} = %.1f, N = %d$ \\' \
              %(exp(logh), Tmax, N) + s, size=23)
ax.set_ylabel(r'Energi per partikkel', size=23)
ax.legend(prop={'size':16})
box = ax.get_position()
ax.set_position([box.x0, box.y0*1.25, box.width, box.height])
plt.show()
