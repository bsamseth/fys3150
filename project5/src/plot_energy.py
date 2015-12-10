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


logh, Tmax, N, check_ =  [float(datafile[:-4].split('_')[i+3]) for i in range(0,4)]

if check_:
    s = 'Her fjerner vi utskytne partikler'
else:
    s = 'Her fjerner vi \_ikke\_ utskytne partikler'


data = loadtxt(datafile)

t_array = [Tmax*i/1000.0 for i in range(len(data[:,0]))]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex = True)
plt.plot(t_array, data[:,0]/N)
plt.title(r'Energi som funksjon av tiden', fontsize='large')
plt.xlabel(r'Tid  [$\tau_{crunch}$] \\ \\' \
           + r'Parametre: $h = %g, T_{max} = %.1f, N = %d$ \\' \
           %(exp(logh), Tmax, N) + s, fontsize='large')
plt.ylabel(r'Energi per partikkel', fontsize='large')
plt.show()

