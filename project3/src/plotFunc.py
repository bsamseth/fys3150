import matplotlib.pyplot as plt
import numpy as np

r = np.linspace(0, 5, 1001)
alpha = 2
psi = np.exp(-alpha * r)


plt.plot(r,psi); 
plt.rc('text', usetex=True)
plt.xlabel(r'$r = r1+r2$', fontsize='large')
plt.ylabel(r'$\Psi$', fontsize='large')
plt.title(r'Plot av $\Psi$ som funksjon av r')

plt.savefig('../fig/psiPlot.png')
plt.show()
