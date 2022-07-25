import numpy as np
import matplotlib.pyplot as plt

n = np.loadtxt("racs.txt")
print(n)
x = np.arange(1, 361)
plt.scatter(x, n)
plt.show()