from re import X
import numpy as np 
import matplotlib.pyplot as plt

def some_fun(x, y):
    return x**2 + 5 * y

x = np.random.rand(10)
y = np.random.rand(10)

# make 2D empty array
z = np.zeros(shape=(len(x), len(y)))

for numx, i in enumerate(x):
    for numy, j in enumerate(y):
        val = some_fun(i, j)
        z[numx,numy] = val

print(z)

#plt.pcolormesh(x, y, z)
#plt.show()
