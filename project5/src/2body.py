import time
from ODESolver import RungeKutta4 as RK4
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Body:
    def __init__(self, mass, initpos, G=6.67408e-11):
        self.G    = G
        self.mass = mass
        self.pos  = initpos
        
    def __call__(self, other):
        return - self.G * self.mass * other.mass / linalg.norm(self.pos - other.pos)**3 * (self.pos - other.pos)


def grav(r1, r2, m1, m2):
    return - G * m1 * m2 / linalg.norm(r1 - r2)**3 * (r1 - r2)
    
m1 = 1. # kg
m2 = 1. # kg
G=1#6.67408e-11

T = 40.
n = 500

h = T/n

t = linspace(0, T, n)
r1 = ndarray((3,n))
r2 = ndarray((3,n))

r1[:,0] = asarray([0,0,0])
r2[:,0] = asarray([5,0,0])

v1 = ndarray((3,n))
v2 = ndarray((3,n))

v1[:,0] = asarray([0,-0.1,0])
v2[:,0] = asarray([0,0.1,0])

dt = t[1] - t[0]
t0 = time.time()

F_prev = grav(r1[:,0], r2[:,0], m1, m2)

for i in xrange(n-1):
    
    r1[:, i+1] = r1[:,i] + v1[:,i] *dt + 0.5*F_prev/m1 * dt**2
    r2[:, i+1] = r2[:,i] + v2[:,i] *dt + 0.5*(-F_prev)/m2 * dt**2

    F = grav(r1[:,i+1], r2[:,i+1], m1, m2)
    
    v1[:, i+1] = v1[:,i] + 0.5* (F_prev + F)/m1 * dt
    v2[:, i+1] = v2[:,i] - 0.5* (F_prev + F)/m2 * dt

    F_prev = F
    
t1 = time.time()
print "time spent per loop = %e" % ((t1-t0) / n)


fig, ax = plt.subplots()
#ax = fig.gca(projection='3d')
plt.plot(r1[0,:], r1[1,:], label='1')
plt.plot(r2[0,:], r2[1,:], label='2')
plt.legend()
# plt.plot(r2[0,:], r2[1,:], label='1')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
#plt.show()


# ========================================
# RK4
# ========================================
def f(u,t):
    v, r = u
    r1,r2 = r

    ny = zeros_like(u)
    ny[0,0] = grav(r1,r2,m1,m2)/m1
    ny[0,1] = -ny[0,0]*m1/m2
    ny[1] = v

    return ny
    
t = linspace(0, T, n)
r1 = ndarray((3,n))
r2 = ndarray((3,n))

r1[:,0] = asarray([0,0,0])
r2[:,0] = asarray([5,0,0])

v1 = ndarray((3,n))
v2 = ndarray((3,n))

v1[:,0] = asarray([0,-0.1,0])
v2[:,0] = asarray([0,0.1,0])

dt = t[1] - t[0]

init_cond = asarray([[v1[:,0],v2[:,0]],[r1[:,0],r2[:,0]]])

print init_cond
t0 = time.time()

problem = RK4(f)
problem.set_initial_condition(init_cond)
u, t = problem.solve(t)


    
    
t1 = time.time()
print "time spent per loop (RK4) = %e" % ((t1-t0) / n)


fig, ax = plt.subplots()
#ax = fig.gca(projection='3d')
plt.plot(r1[0,:], r1[1,:], label='1')
plt.plot(r2[0,:], r2[1,:], label='2')
plt.legend()
# plt.plot(r2[0,:], r2[1,:], label='1')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
plt.show()
