import matplotlib.pyplot as plt
import numpy as np

tau_ = np.linspace(-10,10,1001)
t1 = -tau_ + np.sqrt(1+tau_**2)
t2 = -tau_ - np.sqrt(1+tau_**2)

ones = np.zeros(1001) + np.pi/4
plt.rc('text', usetex=True)
plt.plot(tau_, np.arctan(t1), label=r'$\theta_1 = tan^{-1}(\tau + \sqrt{1+\tau^2}$)')
plt.plot(tau_, np.arctan(t2), label=r'$\theta_2 = tan^{-1}(\tau - \sqrt{1+\tau^2}$)')
plt.plot(tau_,ones, 'r--')
plt.plot(tau_, -ones, 'r--', label=r'$\pm \frac{\pi}{4}$', )
plt.xlabel(r'$ cos 2 \theta / [ \tau ]$', fontsize='large')
plt.ylabel(r'Angle, \theta', fontsize='large')
plt.vlines(0, -2, 2, colors = u'k', linestyles='dashed')# np.pi/2, np.pi/2)
plt.legend()
plt.savefig("fig/theta.png")
plt.show()


