import numpy as np

import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotDensity1d
from pysgpp import createOperationDensityMarginalizeKDE, KernelDensityEstimator, DataMatrix


class NatafTransformation(object):

    def __init__(self, data, sample_type=None, dist=None):
        from pysgpp.extensions.datadriven.uq.dists import Uniform, Beta, SGDEdist, Normal, KDEDist
        from pysgpp.extensions.datadriven.uq.quadrature.marginalization.marginalization import doMarginalize

        # fix stochastic setting
        self.alpha, self.beta = 5., 10.
        self.lwr, self.upr = 0., 1.
        self.normal = Normal(0, 1, -2, 2)
        self.uniform = Uniform(self.lwr, self.upr)
        self.b = Beta(self.alpha, self.beta, self.lwr, self.upr)
        self.dim = data.shape[0]

        if sample_type == 'cbeta':
            # marginalize the density
            opMar = createOperationDensityMargTo1DKDE(dist.dist)
            kdex = KernelDensityEstimator()
            opMar.margToDimX(kdex, 0)
            kdey = KernelDensityEstimator()
            opMar.margToDimX(kdey, 1)

            # set the mean vector and the correlation matrix
            self.x = [KDEDist(kdex.getSamples().array()),
                      KDEDist(kdey.getSamples().array())]
            self.M = np.array([[kdex.mean(), kdey.mean()]]).T
            self.S = dist.corrcoeff()
        else:
            self.x = [self.b, self.b]
            self.M = np.array([[self.b.mean(), self.b.mean()]]).T
            self.S = np.array([[1., 0.],
                               [0., 1.]])

        # compute the correlation matrix from the covariance matrix
        # this is used to transform the results back to the original space
        self.D = np.diag(np.sqrt(np.diag(self.S)))
        # divide the diagonal by the standard deviation of the diagonal elements
        self.D_inverse = np.diag(1. / np.sqrt(np.diag(self.S)))
        self.C = self.D_inverse.dot(self.S.dot(self.D_inverse))

#         fig = plt.figure()
#         plotDensity1d(self.x[0])
#         plotDensity1d(self.b)
#         fig.show()
#
#         fig = plt.figure()
#         plotDensity1d(self.x[1])
#         plotDensity1d(self.b)
#         fig.show()

        # compute cholesky decomposition
        self.L = np.linalg.cholesky(self.C)

        # adjust it according to [Lu ...]
        # nothing needs to be done for uniform <--> uniform
        self.L = self.L
        self.L_inverse = np.linalg.inv(self.L)

        assert abs(np.sum(self.C - self.L.dot(self.L.T))) < 1e-14
        assert abs(np.sum(self.S - self.D.dot(self.L.dot(self.L.T.dot(self.D))))) < 1e-14

    def trans_U_to_X(self, u_vars, x_vars):
        z_vars = np.zeros(u_vars.shape)
        self.trans_U_to_Z(u_vars, z_vars)
        self.trans_Z_to_X(z_vars, x_vars)

    def trans_X_to_U(self, x_vars, u_vars):
        z_vars = np.zeros(u_vars.shape)

        self.trans_X_to_Z(x_vars, z_vars)
        self.trans_Z_to_U(z_vars, u_vars)

    def trans_Z_to_X(self, z_vars, x_vars):
        for i in xrange(self.dim):
            normcdf = self.normal.cdf(z_vars[i])
            scaled_x = self.x[i].ppf(normcdf.reshape(len(normcdf), 1))
            scaled_x = scaled_x.reshape(len(normcdf))
            x_vars[i] = self.lwr + (self.upr - self.lwr) * scaled_x

    def trans_X_to_Z(self, x_vars, z_vars):
        for i in xrange(self.dim):
            betacdf = self.x[i].cdf(x_vars[i].reshape(len(x_vars[i]), 1))
            betacdf = betacdf.reshape(len(betacdf))
            z_vars[i] = self.normal.ppf(betacdf)

    def trans_Z_to_U(self, z_vars, u_vars):
        # decorrelate the variables
        res = self.L_inverse.dot(self.D_inverse.dot(z_vars - self.M))

        # transform to uniform space
        for i, zi in enumerate(res):
            u_vars[i] = self.normal.cdf(zi)

    def trans_U_to_Z(self, u_vars, z_vars):
        # transform to std normal space
        for i, ui in enumerate(u_vars):
            z_vars[i] = self.normal.ppf(ui)

        # apply the correlation
        res = self.D.dot(self.L.dot(z_vars)) + self.M

        # transform to space of correlated normal
        for i, zi in enumerate(res):
            z_vars[i] = zi
