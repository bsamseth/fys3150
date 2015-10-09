import matplotlib.pyplot as plt
import numpy as np

r = np.linspace(0, 5, 1001)
alpha = 2
psi = np.exp(-alpha * r)

plt.plot(r,psi)
plt.show()
