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

logh, Tmax, N, check_ =  [float(datafile[:-4].split('_')[i+3]) for i in range(0,4)]


Total_energy = data[:,0]
Kinetic_energy = data[:,1]
Potential_energy = Total_energy-Kinetic_energy

Equilibrium_index = len(Total_energy)/2
eq = Equilibrium_index


t_array = [Tmax*i/float(len(data[:,0])) for i in range(len(data[:,0]))][eq:]



s = 'Her fjerner vi %sutskytne partikler' % ("\_ikke\_ " if not check_ else "")
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex = True)
fig, ax = plt.subplots()
ax.plot(t_array, -Potential_energy[eq:]/Kinetic_energy[eq:],
        label='Simulation')
ax.plot(t_array, ones_like(Potential_energy[eq:])*2,label='Theorem')
ax.set_xlabel(r'Tid  [$\tau_{crunch}$] \\ \\' \
           + r'Parametre: $h = %g, T_{max} = %.1f, N = %d$ \\' \
           %(exp(logh), Tmax, N) + s, fontsize='large')
ax.set_ylabel(r'$-\langle V\rangle / \langle K\rangle$')
ax.legend()
box = ax.get_position()
ax.set_position([box.x0, box.y0*1.25, box.width, box.height])

ax.annotate(r'$\langle -\frac{\langle V\rangle}{\langle K\rangle}\rangle = %2.3f$'%(sum(-Potential_energy[eq:]/Kinetic_energy[eq:])/eq),
            xy=(1, 0), xycoords='axes fraction', fontsize=16,
            horizontalalignment='right', verticalalignment='bottom')
           
print sum(-Potential_energy[eq:]/Kinetic_energy[eq:])/eq
plt.show()




