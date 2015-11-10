#! /usr/bin/python
"""
Program used to plot P(E) with data produces by MPIisingImproved.cpp
"""

import os
import matplotlib.pyplot as plt
import sys
from time import time
from numpy import linspace, logspace, zeros, where
from math import log10

try:
    n_spins      = int(sys.argv[1])
    T0           = float(sys.argv[2])
    T1           = float(sys.argv[3])
    dT           = float(sys.argv[4])
    MCcycles_N   = int(sys.argv[5])
    randomizer   = int(sys.argv[6])
    noCompute    = False # assume everything needs to be calculated
except:
    print 'Usage: %s n_spins T0 T1 dT MCcycle_N randomizer [--no-compute]' % sys.argv[0]
    sys.exit(1)

if "--no-compute" in sys.argv:
    noCompute = True
    
if not noCompute: os.system("make")  # compile cpp program

t0 = time() # time execution

# set up arrays
T_N   = int(round((T1 - T0) / dT + 1))
T     = linspace(T0, T1, T_N)   
count = zeros(T_N, dtype=int)

# compute/read data
if not noCompute:
    cmd = "mpirun -n 4 MPIisingImproved.x data/out %d %d %g %g %g %d" % (n_spins, MCcycles_N, T0, T1, dT, randomizer)
    print "Running: %s" % cmd
    os.system(cmd)
    
with open("data/out_Energies_%d_%d_%g_%g_%g_%d.dat" % (n_spins, MCcycles_N, T0, T1, dT, randomizer), 'r') as data: # read data
    data = data.readlines()
    E    = zeros((T_N, len(data)/T_N))
    used_counts = [[],[]]
    for j, line in enumerate(data):
        numbers = [float(elm) for elm in line.split()]
        for multiple in range(0, T_N):
            if abs(abs(numbers[1]-T0) - multiple * dT) < 1e-10:
                index = multiple
                break
        else:
            raise IndexError("Could not find a nice multiple for %g in data file on line %d" % (numbers[1], j))    
        # index = int((numbers[1] - T0) // dT)
        E[index, count[index]] = numbers[0]
        count[index] += 1

# done computing
t1 = time()
print 'Complete. Total time spent: %5.3f' % (t1-t0)

# make plots
# making two seperate plots with three quantities per figure
# expectation values in one, and the variances in the other
os.system("mkdir -p ../fig")  #  make sure the dir exists
randomName = "_random" if randomizer else ""
plt.rc('text', usetex=True)
fig, ax = plt.subplots(T_N, sharex=False, sharey=False)
for i in range(T_N):
    ax[i].hist(E[i,:], bins=range(int(min(E[i,:])), int(max(E[i,:]))+1),
            rwidth=1, normed=1, label=r"$T=%g$" % T[i])
    ax[i].legend()
    ax[i].set_ylabel("Prob")
ax[len(ax)-1].set_xlabel("Energy")
plt.show()
