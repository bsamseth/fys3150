#! /usr/bin/python
import numpy as np
import os
import sys
import subprocess
import time

def center(s, width, sep=" "):
    chars = len(s)
    lenspace = (width - chars)/2 
    spaces = sep * lenspace
    return spaces + s + spaces + sep*(chars%2 == 0)

# os.system("make")

N = 100
log_h_list = range(-2,-5,-1)
t_max_list = np.linspace(0.3, 2.5, 10)

width = 6+ 2*len(log_h_list)

print "|   %s%s|" % (center("Verlet", width), center("RK4", width))
print "|   " + "="*(2*width) + "|"
print "|   %s%s|" % (center("step, h", width, sep="_"), center("step, h", width, sep="_"))
print ("|   | " + " ".join([str(h) for h in log_h_list]))*2
# sys.exit(0)
for t_max in t_max_list:
    print "| T=%.2f |" % t_max,
    for solver in ["verlet", "rk4"]:
        for h in [10**i for i in log_h_list]:
            cmd = "./main.x %g %g %g %s" % (h, t_max, N, solver)
            t0 = time.time()
            proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, )
            output=proc.communicate()[0]
            t1 = time.time()
            error = float(output.split()[-1])
            print "(%.1E, %.1f)" % (abs(error), t1-t0),
        print " | ",
    print 
