import okushiri
import numpy as np
import time
import matplotlib.pyplot as plt

# load parameters, can be 2 or 4 dimensional (3?)
x = np.loadtxt('x.txt')
print(type(x))
print(np.shape(x))
# set gridsize, can be up to 256 (or more?) but then is slower
gridsize = np.loadtxt('gridsize.txt')

# reset time
np.savetxt('t.txt', [-1])

# run solver
t = time.time()
y = okushiri.run(x, gridsize)
dt = time.time() - t

# save output
np.savetxt('y.txt', y)
np.savetxt('dt.txt', [dt])

# plt.plot(range(len(y[0])), y[0])
# plt.show()
