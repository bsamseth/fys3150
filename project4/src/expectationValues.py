#! /usr/bin/python
"""
Program used to visualize the expectation values generated
by MPIising.cpp as functions of the number of Monte Carlo cycles.

Usage: ./expectationValues.py n_spins T0 T1 dT MCcycle_min MCcycle_max MCcycle_N randomizer [--no-compute]'

Optional parameter --no-compute can be added as final argument. If so,
all data files are assumed to exits.
"""

import os
import matplotlib.pyplot as plt
import sys
from time import time
from numpy import linspace, logspace, zeros
from math import log10

try:
    n_spins      = int(sys.argv[1])
    T0           = float(sys.argv[2])
    T1           = float(sys.argv[3])
    dT           = float(sys.argv[4])
    MCcycles_min = int(sys.argv[5])
    MCcycles_max = int(sys.argv[6])
    MCcycles_N   = int(sys.argv[7])
    randomizer   = int(sys.argv[8])
    noCompute    = False # assume everything needs to be calculated
except:
    print 'Usage: %s n_spins T0 T1 dT MCcycle_min MCcycle_max MCcycle_N randomizer [--no-compute]' % sys.argv[0]
    sys.exit(1)

if "--no-compute" in sys.argv:
    noCompute = True
    
if not noCompute: os.system("make")  # compile cpp program

t0 = time() # time execution

# set up arrays
T_N      = int(round((T1 - T0) / dT + 1))
T_array  = linspace(T0, T1, T_N)
MCcycles = logspace(log10(MCcycles_min), log10(MCcycles_max), MCcycles_N)
E        = zeros((T_N,MCcycles_N))
Cv       = zeros((T_N,MCcycles_N))
M        = zeros((T_N,MCcycles_N))
Mvar     = zeros((T_N,MCcycles_N))
MvarAbs  = zeros((T_N,MCcycles_N))
Mabs     = zeros((T_N,MCcycles_N))
conf_N   = zeros((T_N,MCcycles_N))

# compute/read data
for i, cycles in enumerate(MCcycles):

    if not noCompute:
        cmd = "mpirun -n 4 MPIising.x data/out %d %d %g %g %g %d" % (n_spins, cycles, T0, T1, dT, randomizer)
        print "Running: %s" % cmd
        os.system(cmd)
    
    with open("data/out_%d_%d_%g_%g_%g_%d.dat" % (n_spins, cycles, T0, T1, dT, randomizer), 'r') as data: # read data
        data = data.readlines()[1:]
        for j, line in enumerate(data):
            numbers = [float(elm) for elm in line.split()[1:]]
            E[j,i], Cv[j,i], M[j,i], Mvar[j,i], MvarAbs[j,i], Mabs[j,i], conf_N[j,i] = numbers

# done computing
t1 = time()
print 'Complete. Total time spent: %5.3f' % (t1-t0)


# make plots
# making two seperate plots with three quantities per figure
# expectation values in one, and the variances in the other
os.system("mkdir -p ../fig")  #  make sure the dir exists
randomName = "_random" if randomizer else ""
fig, ax = plt.subplots(3, sharex=True)
plt.rc('text', usetex=True)
for i in range(T_N):
    ax[0].semilogx(MCcycles, E[i,:], label=r'$T = %g$' % T_array[i])
    ax[1].semilogx(MCcycles, M[i,:], label=r'$T = %g$' % T_array[i])
    ax[2].semilogx(MCcycles, Mabs[i,:], label=r'$T = %g$' % T_array[i])

ax[0].set_ylabel(r"$\langle E \rangle$", size=16)
ax[1].set_ylabel(r"$\langle M \rangle$", size=16)
ax[2].set_ylabel(r"$\langle |M| \rangle$", size=16)
ax[2].set_xlabel(r"Monte Carlo Cycles")
ax[0].set_title(r'Expectation values as functions of Monte Carlo Cycles \\' \
                'with $\text{n_spins T0 T1 dT MCcycle_min MCcycle_max MCcycle_N randomizer} = %s $' \
                %([float(i) for i in sys.argv[1:-1]]))  
ax[len(ax)/2].legend(loc='center left', bbox_to_anchor=(1, 0.5)) # one common legend
plt.savefig("../fig/E_M_Mabs%s.png" % randomName) # might need to change window for all to look good

fig, ax = plt.subplots(3, sharex=True)
for i in range(T_N):
    ax[0].semilogx(MCcycles, Mvar[i,:], label=r'$T = %g$' % T_array[i])
    ax[1].semilogx(MCcycles, MvarAbs[i,:], label=r'$T = %g$' % T_array[i])
    ax[2].semilogx(MCcycles, Cv[i,:], label=r'$T = %g$' % T_array[i])
ax[0].set_ylabel(r"$\sigma_M^2$", size=16)
ax[1].set_ylabel(r"$\sigma_{|M|}^2$", size=16)
ax[2].set_ylabel(r"$\langle C_V \rangle$", size=16)
ax[2].set_xlabel(r"Monte Carlo Cycles")
ax[0].set_title(r'Variances as functions of Monte Carlo Cycles \\' \
                'with $n_spins T0 T1 dT MCcycle_min MCcycle_max MCcycle_N randomizer = %s $' \
                %([float(i) for i in sys.argv[1:-1]])) 
ax[len(ax)/2].legend(loc='center left', bbox_to_anchor=(1, 0.5)) 
plt.savefig("../fig/variances%s.png" % randomName)

fig, ax = plt.subplots()
for i in range(T_N):
    ax.semilogx(MCcycles, conf_N[i,:]/MCcycles, label=r"$T = %g$" % T_array[i])
ax.set_xlabel("Monte Carlo Cycles")
ax.set_ylabel("Configurations")
ax.set_title(r'Number of accepted configurations per Monte Carlo Cycle \\'
             'with $n_spins T0 T1 dT MCcycle_min MCcycle_max MCcycle_N randomizer = %s$' \
             %([float(i) for i in sys.argv[1:-1]])) 
ax.legend()
plt.savefig("../fig/configurations%s.png" % randomName)
plt.show()
