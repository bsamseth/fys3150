#! /usr/bin/python
"""

"""

import os
import matplotlib.pyplot as plt
import sys
from time import time
from numpy import linspace, logspace, zeros
from math import log10

try:
    L            = eval(sys.argv[1])
    T0           = float(sys.argv[2])
    T1           = float(sys.argv[3])
    dT           = float(sys.argv[4])
    MCcycles_N   = int(sys.argv[5])
    randomizer   = int(sys.argv[6])
    noCompute    = False # assume everything needs to be calculated
except:
    print 'Usage: %s [list of L] T0 T1 dT MCcycle_N randomizer [--no-compute]' % sys.argv[0]
    sys.exit(1)

if "--no-compute" in sys.argv:
    noCompute = True
    
if not noCompute: os.system("make")  # compile cpp program

t0 = time() # time execution

# set up arrays
T_N      = int(round((T1 - T0) / dT + 1))
T_array  = linspace(T0, T1, T_N)
L_N      = len(L)
E        = zeros((L_N, T_N))   
Cv       = zeros((L_N, T_N))   
chi      = zeros((L_N, T_N))   
Mabs     = zeros((L_N, T_N))   

# compute/read data
for j, Li in enumerate(L):

    if not noCompute:
        cmd = "mpirun -n 4 MPIising.x data/out %d %d %g %g %g %d" % (Li, MCcycles_N, T0, T1, dT, randomizer)
        print "Running: %s" % cmd
        os.system(cmd)
    
    with open("data/out_%d_%d_%g_%g_%g_%d.dat" % (Li, MCcycles_N, T0, T1, dT, randomizer), 'r') as data: # read data
        data = data.readlines()[1:]
        for i, line in enumerate(data):
            numbers = [float(elm) for elm in line.split()[1:]]
            E[j,i], Cv[j,i], dummy, dummy2, chi[j,i], Mabs[j,i], dummy3 = numbers

# done computing
t1 = time()
print 'Complete. Total time spent: %5.3f' % (t1-t0)


# make plots
# making two seperate plots with three quantities per figure
# expectation values in one, and the variances in the other
os.system("mkdir -p ../fig")  #  make sure the dir exists
randomName = "_random" if randomizer else ""
fig, ax = plt.subplots()
plt.rc('text', usetex=True)
for i in range(L_N):
    ax.plot(T_array, Mabs[i,:], label=r"$L=%g" % L[i])

plt.show()

