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

def center2(s1, s2, width):
    s1 = center(s1,width/2)
    s2 = center(s2,width/2)
    return s1+s2

def center(s, width, sep=""):
    return ('{:%s^%d}' % (sep,width)).format(s)

# os.system("make")

N = 100
log_h_list = range(-2,-5,-1)
t_max_list = np.linspace(0.3, 2.5, 10)

width = "| T=0.00 | " + " ".join(["(0.0E-00, 000.0)" for i in range(len(log_h_list))]) + "  |  " \
        + " ".join(["(0.0E-00, 000.0)" for i in range(len(log_h_list))]) + "  |"
width = len(width)
print "|%s%s|" % (center("Verlet", width/2), center("RK4", width/2 -1))
print "|" + "="*(width-2) + "|"
print "|%s%s|" % (center("log(h)=%s" % ", ".join([str(h) for h in log_h_list]), width/2, sep="_"), center("log(h)=%s" % ", ".join([str(h) for h in log_h_list]), width/2 -1, sep="_"))
# print ("|   | " + )*2


# sys.exit(0)
for t_max in t_max_list:
    print "| T=%.2f |" % t_max,
    for solver in ["verlet", "rk4"]:
        for h in [10**i for i in log_h_list]:
            cmd = "./main.x %g %g %g %d %s" % (h, t_max, N, 0, solver)
            t0 = time.time()
            proc=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, )
            output=proc.communicate()[0]
            t1 = time.time()
            error = float(output.split()[-1])
            print "(%.1E, %05.1f)" % (abs(error), t1-t0),
        print " | ",
    print 
