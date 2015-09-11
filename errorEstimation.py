#! /usr/bin/python
import sys
import numpy as np
from subprocess import call
from scitools.std import *

try:
    N = np.asarray([int(i) for i in sys.argv[1:]])
except IndexError:
    print "Usage: %s list of N" % sys.argv[0]
    sys.exit(1)


eps = np.zeros(len(N))
s = 0
for N_ in N:
    data = open("data/v_solve_N_%d.dat" % N_, 'r')
    data = data.read()
    data = data.split("\n")
    v = zeros(N_)
    for i, line in enumerate(data):
        i -= 1
        if i < 1 or line == "": continue
        v[i] = float(line)
    
    x = linspace(0,1,N_)
    x_exact = linspace(0,1,N_)
    u_exact = x_exact * (exp(-10) - 1) - exp(-10*x_exact)  + 1

    eps[s] = max(np.log10((v-u_exact)/(u_exact+1e-10)))
    s+=1

h_ = 1.0/(N-1)
print h_, eps
semilogx(h_, eps, legend=["$eps$"]); title('Error estimation')
xlabel('$log10(h)$'); ylabel('Max value of $ \epsilon_i$')
raw_input()
