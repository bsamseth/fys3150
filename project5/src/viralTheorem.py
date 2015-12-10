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

Total_energy = data[:,0]
Kinetic_energy = data[:,1]
Potential_energy = Total_energy-Kinetic_energy

Equilibrium_index = len(Total_energy)/2
eq = Equilibrium_index

plt.plot(-Potential_energy[eq:]/Kinetic_energy[eq:])
plt.plot(ones_like(Potential_energy[eq:])*2)
print sum(-Potential_energy[eq:]/Kinetic_energy[eq:])/eq

# plt.plot(2*Kinetic_energy[eq:], label='2K')
# plt.plot(-Potential_energy[eq:], label=' -V')
# plt.xlabel('Tid')
# plt.ylabel('Energy')
# plt.legend()
plt.show()




