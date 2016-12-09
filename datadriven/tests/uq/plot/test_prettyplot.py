#!/usr/bin/env python
"""
pcolormesh uses a QuadMesh, a faster generalization of pcolor, but
with some restrictions.

This demo illustrates a bug in quadmesh with masked data.
"""

from matplotlib import cm, colors
from numpy import ma

import numpy as np
import matplotlib.pyplot as plt


n = 12
x = np.linspace(-1.5, 1.5, n)
y = np.linspace(-1.5, 1.5, n * 2)
X, Y = np.meshgrid(x, y)
Qx = np.cos(Y) - np.cos(X)
Qz = np.sin(Y) + np.sin(X)
Qx = (Qx + 1.1)
Z = np.sqrt(X ** 2 + Y ** 2) / 5
Z = (Z - Z.min()) / (Z.max() - Z.min())

# The color array can include masked values:
# Zm = ma.masked_where(np.abs(Qz) < 0.5 * np.amax(Qz), Z)


fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_axis_bgcolor("#bdb76b")
# ax.pcolormesh(Qx, Qz, Z, shading='gouraud')
# ax.set_title('Without masked values')
im = plt.imshow(Z, interpolation='bilinear', origin='lower',
                cmap=cm.coolwarm, extent=(-3, 3, -2, 2))
levels = np.arange(-1.2, 1.6, 0.2)
CS = plt.contour(Z, levels,
                 origin='lower',
                 linewidths=2,
                 extent=(-3, 3, -2, 2))

plt.contour(Qx, Qz, Z, 6, colors='k')

plt.show()
