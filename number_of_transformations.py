#! /usr/bin/python
import subprocess
import sys
from numpy import *
import matplotlib.pyplot as plt

subprocess.call(["make", "solve_lower_3_states"])  # compile cpp program

try:
    rho_max = float(sys.argv[1])
    N_list = [int(arg) for arg in sys.argv[1:]]
except:
    print "Usage: %s rho_max N1 N2 N3 ..."
    sys.exit(0)

number_of_transformations_list = []
for N in N_list:
    mainProcess = subprocess.Popen(['./solve_lower_3_states.x', '%g' % N , '%g' % rho_max], stdout=subprocess.PIPE)
    out, err = mainProcess.communicate()
    number_of_transformations = int(out.split()[-3])
    number_of_transformations_list.append(number_of_transformations)


fitted_curve_coeff = polyfit(N_list, number_of_transformations_list, 2)
x = linspace(N_list[0], N_list[-1], 501)
fitted_curve = fitted_curve_coeff[0] * x**2 + fitted_curve_coeff[1] * x + fitted_curve_coeff[2]

print fitted_curve_coeff

x_optimal = linspace(N_list[0], N_list[-1], 501)
y_optimal = 0.5 * x_optimal * ( x_optimal-1)

fig, ax = plt.subplots()
ax.hold(True)
ax.plot(N_list, number_of_transformations_list, label='Actual')
ax.plot(x, fitted_curve, 'r-', label='Trend curve')
ax.plot(x_optimal, y_optimal, label='Optimal curve')
ax.set_title(r"Number of similarity transformations as funciton of $n_{step}$")
ax.set_xlabel(r"$n_{step}$", size=16)
ax.set_ylabel(r"Number of similarity transformations")
ax.legend(loc='upper left')
plt.savefig("fig/number_of_transformations_trend_plot.jpg")
plt.show()
