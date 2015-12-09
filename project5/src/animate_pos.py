import sys, os
import numpy as np
#from scipy import integrate

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation



try:
    datafile = sys.argv[1]
except Exception:
    print "Usage: %s datafile" % sys.argv[0]
    sys.exit(1)

data = np.loadtxt(datafile)

n_bodies = len(data[0])/3
n_timesteps = len(data)

r = np.empty((n_bodies,n_timesteps,3))

print data.shape, r.shape
for i in range(n_bodies):
    r[i] = data[:,i*3:i*3+3]




# Set up figure & 3D axis for animation
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
# ax.axis('off')

# choose a different color for each trajectory
colors = plt.cm.jet(np.linspace(0, 1, n_bodies))

# set up lines and points
lines = sum([ax.plot([], [], [], '-', c=c)
             for c in colors], [])
pts = sum([ax.plot([], [], [], 'o', c=c)
           for c in colors], [])

# prepare the axes limits
xmin_max = (np.min(r[:,:,0])/2, np.max(r[:,:,0])/2)
ymin_max = (np.min(r[:,:,1])/2, np.max(r[:,:,1])/2)
zmin_max = (np.min(r[:,:,2])/2, np.max(r[:,:,2])/2)
ax.set_xlim(xmin_max)
ax.set_ylim(ymin_max)
ax.set_zlim(zmin_max)

# set point-of-view: specified by (altitude degrees, azimuth degrees)
ax.view_init(30, 0)

# initialization function: plot the background of each frame
def init():
    for line, pt in zip(lines, pts):
        line.set_data([], [])
        line.set_3d_properties([])

        pt.set_data([], [])
        pt.set_3d_properties([])
    return lines + pts

# animation function.  This will be called sequentially with the frame number
def animate(i):
    # we'll step two time-steps per frame.  This leads to nice results.
    i = (5 * i) % r.shape[1]

    for line, pt, ri in zip(lines, pts, r):
        # x, y, z = ri[:i].T
        x, y, z = ri[i-1].T
        line.set_data(x, y)
        line.set_3d_properties(z)

        pt.set_data(x, y)
        # pt.set_data(x[-1:], y[-1:])
        pt.set_3d_properties(z)
        # pt.set_3d_properties(z[-1:])

    ax.legend(['t = %g' % (i/float(n_timesteps))])
    #ax.view_init(30, 0.01 *0.3 * i )
    fig.canvas.draw()
    return lines + pts

# instantiate the animator.
#mywriter = animation.FFMpegWriter()
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=n_timesteps, interval=10, blit=True)

# Save as mp4. This requires mplayer or ffmpeg to be installed
anim.save('solarsystem.mp4', fps=15)#, extra_args=['-vcodec', 'libx264'])

plt.show()
