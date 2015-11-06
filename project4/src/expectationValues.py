#! /usr/bin/python
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
call(["make"])  # compile cpp program

fig = plt.figure()
ax = fig.gca()
plt.rc("text", usetex=True)
ax.hold('on')

N = 3
average_E = np.zeros((N,2))
average_Cv = np.zeros((N,2))
average_absM = np.zeros((N,2))
average_M = np.zeros((N,2))
average_sus = np.zeros((N,2))
average_abssus = np.zeros((N,2))
cycles = np.zeros(N)

i = 0
while i < N: 
    cycles[i] = 10000*(i+1)
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
        line = line.split()
        average_E[i, j] = float(line[1])
        average_Cv[i, j] = float(line[2])
        average_M[i, j] = float(line[3])
        average_sus[i, j] = float(line[4])
        average_abssus[i, j] = float(line[5])
        average_absM[i, j] = float(line[6])
    i+=1

plt.rc('text', usetex=True)
plt.subplot(211)
plt.plot(cycles, average_E[:,0])
plt.plot(cycles, average_E[:,1])
plt.ylabel('')
plt.title('Energy')
plt.legend(['T=1', 'T=2.4'])

plt.subplot(212)
plt.plot(cycles, average_absM[:,0])
plt.plot(cycles, average_absM[:,1])
plt.xlabel('Number of MC cycles'); 
plt.ylabel('')
plt.title(r'$|M|$')
plt.legend(['T=1', 'T=2.4'])
plt.show()
