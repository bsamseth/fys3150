#! /usr/bin/python
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
call(["make"])  # compile cpp program

fig = plt.figure()
ax = fig.gca()
plt.rc("text", usetex=True)
ax.hold('on')

N = 1
average_E = np.zeros(N)
average_Cv = np.zeros(N)
average_absM = np.zeros(N)
average_M = np.zeros(N)
average_sus = np.zeros(N)
cycles = np.zeros(N)

i = 0
while i < N: 
    cycles[i] = 1000*(i+1)
    nspins = 20
    T0 = 1.0
    T1 = 2.4
    dT = 1.4

    call(["mpirun", "-n", "4", "MPIising.x", "data/out",
          "%d" %nspins,  "%d" %cycles[i],  "%g" %T0,  "%g" %T1,  "%g" %dT])
    data = open("data/out_%d_%d_%g_%g_%g.dat" %(nspins, cycles[i], T0, T1, dT), 'r')  # read data back in
    data = data.read()
    data = data.split("\n")
    
    for j, line in enumerate(data[:-2]):
        line.split()
        #average_E[i] = float(line[1])
        print line[1]
        print j, line
    i+=1
'''plt.plot(cycles, average_E)
plt.show()'''
