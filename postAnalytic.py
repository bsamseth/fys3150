#! /usr/bin/python
import sys
import numpy as np
from subprocess import call
from scitools.std import *

try:
    N = int(sys.argv[1])
except IndexError:
    print "Usage: %s dimensionN" % sys.argv[0]
    sys.exit(1)

call(["make"])
call(["./main.x", "%d" % N])
    
data = open("data/v_solve_N_%d.dat" % N, 'r')
data = data.read()
data = data.split("\n")
N = int(data[0])
v = zeros(N)
for i, line in enumerate(data):
    i -= 1
    if i < 1 or line == "": continue
    v[i] = float(line)

x = linspace(0,1,N)
x_exact = linspace(0,1,101)
u_exact = x_exact * (exp(-10) - 1) - exp(-10*x_exact)  + 1


plot(x, v, x_exact, u_exact, legend=["$v_i$", "$u(x)$"])
raw_input()


