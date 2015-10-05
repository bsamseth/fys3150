#! /usr/bin/python
import subprocess
import sys
from numpy import *
import matplotlib.pyplot as plt

try:
    nstep_list = []
    rho_max_list = []
    omega_list = []
    if len(sys.argv) - 1 < 3:
        raise IndexError
    for i in range(1, (len(sys.argv)-1), 3):
        nstep_list.append( int(sys.argv[i]) )
        rho_max_list.append( float(sys.argv[i+1]) )
        omega_list.append( float(sys.argv[i+2]) )
except IndexError:
    print "Usage: %s nstep1 rho_max1 omega1 nstep2 ..."
    sys.exit(1)
except TypeError:
    print "Usage: %s nstep1 rho_max1 omega1 nstep2 ..."
    print "All arguments shoudl be numbers, and nstep must be an integer."
    sys.exit(1)

subprocess.call(["make", "interactingElectrons"])  # compile cpp program to exec


for interacting in [0,1]:
    # make figure to put plot in
    number_of_plots = len(omega_list)
    fig, ax = plt.subplots(nrows=number_of_plots)
    plt.hold(True)
    for k in xrange(len(omega_list)):
        omega = omega_list[k] # for shorthand
        nstep = nstep_list[k]
        rho_max = rho_max_list[k]
        
        subprocess.call(["./interactingElectrons.x", str(interacting), str(nstep), str(rho_max), str(omega)])
        noninteractiong = "Non" if interacting == 0 else ""
        dataFile = open("data/%sinteractingElectrons_omegar_%g_nstep_%d_rhomax_%g.dat" % (noninteractiong, omega, nstep, rho_max), 'r')
        data = dataFile.read()
        data = data.split('\n')
        eigenstates = []
        for i, line in enumerate(data):
            if i < 3:
                continue # skip eigenvaules
            if line.strip() == "":
                continue # if line is blank
            # each line is the eigenstate elements, not including the fixed endpoints
            # add 0 as element before and after comuted state.
            elements = [0] + [float(word) for word in line.split()] + [0] 
            eigenstates.append(asarray(elements))
            
        # plot the wavefunction
        h = rho_max / float(nstep)
        rho = asarray([ i*h for i in xrange(nstep+1)])

        cax = ax if len(omega_list) < 2 else ax[k]
        if k == 0:
            cax.set_title(r"Three lowest states, %sinteracting" % (noninteractiong.lower()))
        for n, state in enumerate(eigenstates):
            cax.plot(rho, abs(state)**2, label=r"$\omega_r=%g, n=%d$" % (omega,n))
            cax.set_xlabel(r"$\rho$", size=16)
            # cax.set_ylabel(r"$|\psi|^2$", size=16)
        # Shrink current axis by 20%
        box = cax.get_position()
        cax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        cax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
