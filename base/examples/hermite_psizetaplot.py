from matplotlib import pyplot as plt
import numpy as np
import pysgpp


def testplot():

    X = np.linspace(0, 1, 900)

    psiBase = pysgpp.SPsiHermiteBase()
    zetaBase = pysgpp.SZetaHermiteBase()

    psi_values = [psiBase.eval(1, 1, x) for x in X]
    zeta_values = [zetaBase.eval(1, 1, x) for x in X]

    return X, psi_values, zeta_values


def psi(x):
    if (x < 0):
        return -2 * x ** 3 - 3 * x ** 2 + 1
    else:
        return 2 * x ** 3 - 3 * x ** 2 + 1


def zeta(x):
    if(x < 0):
        return x ** 3 + 2 * x ** 2 + x
    if(x >= 0):
        return x ** 3 - 2 * x ** 2 + x


def uniformZetaPsi():
    X = np.linspace(-1, 1, 900)
    psi_values = [psi(x) for x in X]
    zeta_values = [zeta(x) for x in X]

    return X, psi_values, zeta_values


X, psi_values, zeta_values = testplot()


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(X,  psi_values, label='$\psi$')
ax.plot(X, zeta_values, label="$\zeta$")
plt.legend()

# x axis "labels"
ax.set_xticks([0, 0.5, 1])
ax.set_xticklabels([0, 1, 2])
ax.set_xlabel("index")


# copy paste for the 0 centered axis
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.show()
