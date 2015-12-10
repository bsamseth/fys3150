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
plt.plot(data[:,0])
plt.show()

