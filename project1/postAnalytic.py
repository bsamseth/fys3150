#! /usr/bin/python
from subprocess import call
from numpy import *
import matplotlib.pyplot as plt
call(["make"])  # compile cpp program

fig = plt.figure()
ax = fig.gca()
plt.rc("text", usetex=True)
ax.hold('on')
for N in [10, 100, 1000]:  # loop over the values of N we want to test

    call(["./main.x", "%d" % N])  # compute and write data

    data = open("data/v_solve_N_%d.dat" % N, 'r')  # read data back in
    data = data.read()
    data = data.split("\n")
    N = int(data[0])
    v = zeros(N)
    for i, line in enumerate(data):
        i -= 1
        if i < 1 or line == "": continue
        v[i] = float(line)

    # now plot the computed data
    x = linspace(0, 1, N)
    ax.plot(x, v, label="$v_i (N=%d)$" % N)

# add closed-form solution to plot for reference
x_exact = linspace(0, 1, 101)
u_exact = x_exact * (exp(-10) - 1) - exp(-10*x_exact) + 1
ax.plot(x_exact, u_exact, label="$u(x)$")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$u(x)$")
ax.legend()
ax.set_title("Solution of the 1-D Poisson's equation")
plt.savefig("fig/b.png")
plt.show()
