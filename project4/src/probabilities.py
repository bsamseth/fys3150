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
count = zeros(T_N)
E     = zeros((T_N, MCcycles_N)) + 123


# compute/read data
if not noCompute:
    cmd = "mpirun -n 4 MPIisingImproved.x data/out %d %d %g %g %g %d" % (n_spins, MCcycles_N, T0, T1, dT, randomizer)
    print "Running: %s" % cmd
    os.system(cmd)
    
with open("data/out_Energies_%d_%d_%g_%g_%g_%d.dat" % (n_spins, MCcycles_N, T0, T1, dT, randomizer), 'r') as data: # read data
    data = data.readlines()
    for j, line in enumerate(data):
        numbers = [float(elm) for elm in line.split()]
        index = int((numbers[1] - T0) // dT)
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
fig, ax = plt.subplots()
for i in range(T_N):
    ax.hist(E[i,:], bins=50, normed=1, label=r"$T=%g$" % T[i], alpha=1 - 0.5*i)
ax.legend()
ax.set_xlabel("Energy")
ax.set_ylabel("Prob")
plt.show()
