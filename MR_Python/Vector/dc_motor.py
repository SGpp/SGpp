import sympy
from scipy.integrate import odeint
import numpy as np
import pysgpp

    # the academic dc_motor model we use in cooperation with Bosch
    # Input:
    # R: resistance (default constant value = 9.0)
    # L: inductance (default constant value = 0.11)
    # cm: motor constant 1 (default constant value = 0.5)
    # cg: motor constant 2 (default constant value = 3.0)
    # J: torque of inertia (default constant value = 0.1)
    # d: friction (default constant value = 0.1)
    # Output: current I/x_1 and torque W/x_2
    #    L * \dot{x}_1 &= -R x_1 - c_m * x_2 + 12
    #    J * \dot{x}_2 &= c_g * x_1 - d * x_2.

class dc_motor_analytical_I():
    # Analytical solution for the ODE system, from the Master thesis of Livia Stohrer
    # returns the curent I
    def __init__(self, dim, numTimeSteps):
        self.dim = dim
        self.numTimeSteps = numTimeSteps
        self.pdfs = pysgpp.DistributionsVector()
        # parameter order 'sobolfive' cg, R, cm, L, J, d
        # -> parameter intervals are true values +/- 5%
        cg_true = 3.0
        R_true = 9.0
        cm_true = 0.5
        L_true = 0.11
        J_true = 0.1
        d_true = 0.1
        paramDeviation = 0.05
        lDev = 1 - paramDeviation
        rDev = 1 + paramDeviation
        self.pdfs.push_back(pysgpp.DistributionUniform(
            lDev * cg_true, rDev * cg_true))
        self.pdfs.push_back(pysgpp.DistributionUniform(
            lDev * R_true, rDev * R_true))
        if dim > 2:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * cm_true, rDev * cm_true))
        if dim > 3:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * L_true, rDev * L_true))
        if dim > 4:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * J_true, rDev * J_true))
        if dim > 5:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * d_true, rDev * d_true))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub

    def getName(self):
        return "dc_motor_analytical_I{}D_{}T".format(self.getDim(), self.getOut())

    def getDim(self):
        return self.dim

    def getOut(self):
        return self.numTimeSteps

    def eval(self, x):
        # true parameters are overwritten, depending on the desired dimensionality
        cm = 0.5
        L = 0.11
        J = 0.1
        d = 0.1

        cg = x[0]
        R = x[1]
        if self.dim > 2:
            cm = x[2]
        if self.dim > 3:
            L = x[3]
        if self.dim > 4:
            J = x[4]
        if self.dim > 5:
            d = x[5]

        tStart = 0.0
        tEnd = 1.0
        dt = ((tEnd - tStart) / (self.numTimeSteps-1))
        timeSteps = np.arange(tStart, tEnd+dt, dt)

        V = 12
        mu = np.sqrt((J*R-L*d)**2-4*cm*cg*J*L)
        lambda1 = -(L*d-mu+J*R)/(2*J*L)
        lambda2 = -(L*d+mu+J*R)/(2*J*L)
        v1 = [d/cg-(L*d-mu+J*R)/(2*L*cg), 1]
        v2 = [d/cg-(L*d+mu+J*R)/(2*L*cg), 1]
        c1 = -V/(L*lambda1*(v2[0]-v1[0]))
        c2 = V/(L*lambda2*(v2[0]-v1[0]))

        I = np.zeros(self.numTimeSteps)
        for i, t in enumerate(timeSteps):
            I[i] = c1 * np.exp(lambda1*t)*v1[0]+c2*np.exp(lambda2*t)*v2[0] + \
                (V*(lambda1*v2[0]-lambda2*v1[0])) / \
                (L*lambda1*lambda2*(v1[0]-v2[0]))
        return I

    def evalJacobian(self, x):
        # define gradients
        R = sympy.Symbol('R', real=True)
        L = sympy.Symbol('L', real=True)
        cm = sympy.Symbol('cm', real=True)
        cg = sympy.Symbol('cg', real=True)
        J = sympy.Symbol('J', real=True)
        d = sympy.Symbol('d', real=True)

        t = sympy.Symbol('t', real=True)
        V = 12

        mu = sympy.sqrt((J*R-L*d)**2-4*cm*cg*J*L)

        lambda1 = -(L*d-mu+J*R)/(2*J*L)
        lambda2 = -(L*d+mu+J*R)/(2*J*L)
        v11 = d/cg - (L*d-mu+J*R)/(2*L*cg)
        v21 = d/cg - (L*d+mu+J*R)/(2*L*cg)
        c1 = -V/(L*lambda1*(v21-v11))
        c2 = V/(L*lambda2*(v21-v11))

        I = c1 * sympy.exp(lambda1*t)*v11+c2*sympy.exp(lambda2*t)*v21 + \
            (V*(lambda1*v21-lambda2*v11))/(L*lambda1*lambda2*(v11-v21))

        dIdR = I.diff(R)
        dIdL = I.diff(L)
        dIdcm = I.diff(cm)
        dIdcg = I.diff(cg)
        dIdJ = I.diff(J)
        dIdd = I.diff(d)

        # evaluate
        tStart = 0.0
        tEnd = 1.0
        dt = ((tEnd - tStart) / (self.numTimeSteps-1))
        timeSteps = np.arange(tStart, tEnd+dt, dt)

        # true parameters are overwritten, depending on the desired dimensionality
        cm_eval = 0.5
        L_eval = 0.11
        J_eval = 0.1
        d_eval = 0.1

        cg_eval = x[0]
        R_eval = x[1]
        if self.dim > 2:
            cm_eval = x[2]
        if self.dim > 3:
            L_eval = x[3]
        if self.dim > 4:
            J_eval = x[4]
        if self.dim > 5:
            d_eval = x[5]

        jacobian = np.zeros((self.numTimeSteps, self.getDim()))
        for j, time in enumerate(timeSteps):
            dIdcg_eval = dIdcg.subs(
                {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
            jacobian[j, 0] = float(sympy.N(dIdcg_eval))
            dIdR_eval = dIdR.subs(
                {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
            jacobian[j, 1] = float(sympy.N(dIdR_eval))
            if self.dim > 2:
                dIdcm_eval = dIdcm.subs(
                    {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
                jacobian[j, 2] = float(sympy.N(dIdcm_eval))
            if self.dim > 3:
                dIdL_eval = dIdL.subs(
                    {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
                jacobian[j, 3] = float(sympy.N(dIdL_eval))
            if self.dim > 4:
                dIdJ_eval = dIdJ.subs(
                    {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
                jacobian[j, 4] = float(sympy.N(dIdJ_eval))
            if self.dim > 5:
                dIdd_eval = dIdd.subs(
                    {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
                jacobian[j, 5] = float(sympy.N(dIdd_eval))

        return jacobian

    def getDistributions(self):
        return self.pdfs

    def getMeans(self):
        print("Means unknown")
        means = np.zeros(self.getOut())
        return means

    def getVars(self):
        print("Vars unknown")
        variances = np.zeros(self.getOut())
        return variances

    def getIntegrals(self):
        print("Integrals unknown")
        integrals = np.zeros(self.getOut())
        return integrals


class dc_motor_analytical_W():
    # Analytical solution for the ODE system, from the Master thesis of Livia Stohrer
    # returns the torque w
    def __init__(self, dim, numTimeSteps):
        self.dim = dim
        self.numTimeSteps = numTimeSteps
        self.pdfs = pysgpp.DistributionsVector()
        # parameter order 'sobolfive' cg, R, cm, L, J, d
        # -> parameter intervals are true values +/- 5%
        cg_true = 3.0
        R_true = 9.0
        cm_true = 0.5
        L_true = 0.11
        J_true = 0.1
        d_true = 0.1
        paramDeviation = 0.05
        lDev = 1 - paramDeviation
        rDev = 1 + paramDeviation
        self.pdfs.push_back(pysgpp.DistributionUniform(
            lDev * cg_true, rDev * cg_true))
        self.pdfs.push_back(pysgpp.DistributionUniform(
            lDev * R_true, rDev * R_true))
        if dim > 2:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * cm_true, rDev * cm_true))
        if dim > 3:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * L_true, rDev * L_true))
        if dim > 4:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * J_true, rDev * J_true))
        if dim > 5:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * d_true, rDev * d_true))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub

    def getName(self):
        return "dc_motor_analytical_W{}D_{}T".format(self.getDim(), self.getOut())

    def getDim(self):
        return self.dim

    def getOut(self):
        return self.numTimeSteps

    # Analytical solution for the ODE system, from the Master thesis of Livia Stohrer
    def eval(self, x):
        # true parameters are overwritten, depending on the desired dimensionality
        cm = 0.5
        L = 0.11
        J = 0.1
        d = 0.1

        cg = x[0]
        R = x[1]
        if self.dim > 2:
            cm = x[2]
        if self.dim > 3:
            L = x[3]
        if self.dim > 4:
            J = x[4]
        if self.dim > 5:
            d = x[5]

        tStart = 0.0
        tEnd = 1.0
        dt = ((tEnd - tStart) / (self.numTimeSteps-1))
        timeSteps = np.arange(tStart, tEnd+dt, dt)

        V = 12
        mu = np.sqrt((J*R-L*d)**2-4*cm*cg*J*L)
        lambda1 = -(L*d-mu+J*R)/(2*J*L)
        lambda2 = -(L*d+mu+J*R)/(2*J*L)
        v1 = [d/cg-(L*d-mu+J*R)/(2*L*cg), 1]
        v2 = [d/cg-(L*d+mu+J*R)/(2*L*cg), 1]
        c1 = -V/(L*lambda1*(v2[0]-v1[0]))
        c2 = V/(L*lambda2*(v2[0]-v1[0]))

        W = np.zeros(self.numTimeSteps)
        for i, t in enumerate(timeSteps):
            W[i] = c1*np.exp(lambda1*t) + c2*np.exp(lambda2*t) + \
                (V*(lambda1-lambda2))/(L*lambda1*lambda2*(v1[0]-v2[0]))
        return W

    def evalJacobian(self, x):
        # define gradients
        R = sympy.Symbol('R', real=True)
        L = sympy.Symbol('L', real=True)
        cm = sympy.Symbol('cm', real=True)
        cg = sympy.Symbol('cg', real=True)
        J = sympy.Symbol('J', real=True)
        d = sympy.Symbol('d', real=True)

        t = sympy.Symbol('t', real=True)
        V = 12

        mu = sympy.sqrt((J*R-L*d)**2-4*cm*cg*J*L)

        lambda1 = -(L*d-mu+J*R)/(2*J*L)
        lambda2 = -(L*d+mu+J*R)/(2*J*L)
        v11 = d/cg - (L*d-mu+J*R)/(2*L*cg)
        v21 = d/cg - (L*d+mu+J*R)/(2*L*cg)
        c1 = -V/(L*lambda1*(v21-v11))
        c2 = V/(L*lambda2*(v21-v11))

        W = c1*sympy.exp(lambda1*t)+c2*sympy.exp(lambda2*t) + \
            (V*(lambda1-lambda2))/(L*lambda1*lambda2*(v11-v21))

        dWdR = W.diff(R)
        dWdL = W.diff(L)
        dWdcm = W.diff(cm)
        dWdcg = W.diff(cg)
        dWdJ = W.diff(J)
        dWdd = W.diff(d)

        # evaluate
        tStart = 0.0
        tEnd = 1.0
        dt = ((tEnd - tStart) / (self.numTimeSteps-1))
        timeSteps = np.arange(tStart, tEnd+dt, dt)

        # true parameters are overwritten, depending on the desired dimensionality
        cm_eval = 0.5
        L_eval = 0.11
        J_eval = 0.1
        d_eval = 0.1

        cg_eval = x[0]
        R_eval = x[1]
        if self.dim > 2:
            cm_eval = x[2]
        if self.dim > 3:
            L_eval = x[3]
        if self.dim > 4:
            J_eval = x[4]
        if self.dim > 5:
            d_eval = x[5]

        jacobian = np.zeros((self.numTimeSteps, self.getDim()))
        for j, time in enumerate(timeSteps):
            dWdcg_eval = dWdcg.subs(
                {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
            jacobian[j, 0] = float(sympy.N(dWdcg_eval))
            dWdR_eval = dWdR.subs(
                {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
            jacobian[j, 1] = float(sympy.N(dWdR_eval))
            if self.dim > 2:
                dWdcm_eval = dWdcm.subs(
                    {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
                jacobian[j, 2] = float(sympy.N(dWdcm_eval))
            if self.dim > 3:
                dWdL_eval = dWdL.subs(
                    {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
                jacobian[j, 3] = float(sympy.N(dWdL_eval))
            if self.dim > 4:
                dWdJ_eval = dWdJ.subs(
                    {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
                jacobian[j, 4] = float(sympy.N(dWdJ_eval))
            if self.dim > 5:
                dWdd_eval = dWdd.subs(
                    {t: time, cg: cg_eval, R: R_eval, cm: cm_eval, L: L_eval, J: J_eval, d: d_eval})
                jacobian[j, 5] = float(sympy.N(dWdd_eval))

        return jacobian

    def getDistributions(self):
        return self.pdfs

    def getMeans(self):
        print("Means unknown")
        means = np.zeros(self.getOut())
        return means

    def getVars(self):
        print("Vars unknown")
        variances = np.zeros(self.getOut())
        return variances

    def getIntegrals(self):
        print("Integrals unknown")
        integrals = np.zeros(self.getOut())
        return integrals

def dc_motor_system(state, t, R=9.0, L=0.11, cm=0.5, cg=3.0, J=0.1, d=0.1):
    x1, x2 = state
    V = 12  # usually voltage of 12 or 24. Could be made a 7th parameter
    d_x1 = (-R * x1 - cm * x2 + V) / L
    d_x2 = (cg * x1 - d * x2) / J
    return[d_x1, d_x2]

def dc_motor_no_load(timeSteps, R=9.0, L=0.11, cm=0.5, cg=3.0, J=0.1, d=0.1):
    # # initial values for current x1 and torque x2
    init_state = [0, 0]
    # # solve ode system
    sol = odeint(dc_motor_system, init_state,
                 timeSteps, args=(R, L, cm, cg, J, d))
    return sol

class dc_motor_ode_I():
    # Numerically solve the ODE system and use this approximate solution
    # returns the curent I
    def __init__(self, dim, numTimeSteps):
        self.dim = dim
        self.numTimeSteps = numTimeSteps
        self.pdfs = pysgpp.DistributionsVector()
        # parameter order 'sobolfive' cg, R, cm, L, J, d
        # -> parameter intervals are true values +/- 5%
        cg_true = 3.0
        R_true = 9.0
        cm_true = 0.5
        L_true = 0.11
        J_true = 0.1
        d_true = 0.1
        paramDeviation = 0.05
        lDev = 1 - paramDeviation
        rDev = 1 + paramDeviation
        self.pdfs.push_back(pysgpp.DistributionUniform(
            lDev * cg_true, rDev * cg_true))
        self.pdfs.push_back(pysgpp.DistributionUniform(
            lDev * R_true, rDev * R_true))
        if dim > 2:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * cm_true, rDev * cm_true))
        if dim > 3:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * L_true, rDev * L_true))
        if dim > 4:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * J_true, rDev * J_true))
        if dim > 5:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * d_true, rDev * d_true))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub

    def getName(self):
        return "dc_motor_ode_I{}D_{}T".format(self.getDim(), self.getOut())

    def getDim(self):
        return self.dim

    def getOut(self):
        return self.numTimeSteps

    def eval(self, x):
        # true parameters are overwritten, depending on the desired dimensionality
        cm = 0.5
        L = 0.11
        J = 0.1
        d = 0.1

        cg = x[0]
        R = x[1]
        if self.dim > 2:
            cm = x[2]
        if self.dim > 3:
            L = x[3]
        if self.dim > 4:
            J = x[4]
        if self.dim > 5:
            d = x[5]

        tStart = 0.0
        tEnd = 1.0
        dt = ((tEnd - tStart) / (self.numTimeSteps-1))
        timeSteps = np.arange(tStart, tEnd+dt, dt)

        sol = dc_motor_no_load(timeSteps, R, L, cm, cg, J, d)
        I = sol[:,0]
        return I

    def evalJacobian(self, x):
        print("Jacobian unknown")
        return 0

    def getDistributions(self):
        return self.pdfs

    def getMeans(self):
        print("Means unknown")
        means = np.zeros(self.getOut())
        return means

    def getVars(self):
        print("Vars unknown")
        variances = np.zeros(self.getOut())
        return variances

    def getIntegrals(self):
        print("Integrals unknown")
        integrals = np.zeros(self.getOut())
        return integrals


class dc_motor_ode_W():
    # Numerically solve the ODE system and use this approximate solution
    # returns the curent I
    def __init__(self, dim, numTimeSteps):
        self.dim = dim
        self.numTimeSteps = numTimeSteps
        self.pdfs = pysgpp.DistributionsVector()
        # parameter order 'sobolfive' cg, R, cm, L, J, d
        # -> parameter intervals are true values +/- 5%
        cg_true = 3.0
        R_true = 9.0
        cm_true = 0.5
        L_true = 0.11
        J_true = 0.1
        d_true = 0.1
        paramDeviation = 0.05
        lDev = 1 - paramDeviation
        rDev = 1 + paramDeviation
        self.pdfs.push_back(pysgpp.DistributionUniform(
            lDev * cg_true, rDev * cg_true))
        self.pdfs.push_back(pysgpp.DistributionUniform(
            lDev * R_true, rDev * R_true))
        if dim > 2:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * cm_true, rDev * cm_true))
        if dim > 3:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * L_true, rDev * L_true))
        if dim > 4:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * J_true, rDev * J_true))
        if dim > 5:
            self.pdfs.push_back(pysgpp.DistributionUniform(
                lDev * d_true, rDev * d_true))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub

    def getName(self):
        return "dc_motor_ode_W{}D_{}T".format(self.getDim(), self.getOut())

    def getDim(self):
        return self.dim

    def getOut(self):
        return self.numTimeSteps

    def eval(self, x):
        # true parameters are overwritten, depending on the desired dimensionality
        cm = 0.5
        L = 0.11
        J = 0.1
        d = 0.1

        cg = x[0]
        R = x[1]
        if self.dim > 2:
            cm = x[2]
        if self.dim > 3:
            L = x[3]
        if self.dim > 4:
            J = x[4]
        if self.dim > 5:
            d = x[5]

        tStart = 0.0
        tEnd = 1.0
        dt = ((tEnd - tStart) / (self.numTimeSteps-1))
        timeSteps = np.arange(tStart, tEnd+dt, dt)

        sol = dc_motor_no_load(timeSteps, R, L, cm, cg, J, d)
        W = sol[:, 1]
        return W

    def evalJacobian(self, x):
        print("Jacobian unknown")
        return 0

    def getDistributions(self):
        return self.pdfs

    def getMeans(self):
        print("Means unknown")
        means = np.zeros(self.getOut())
        return means

    def getVars(self):
        print("Vars unknown")
        variances = np.zeros(self.getOut())
        return variances

    def getIntegrals(self):
        print("Integrals unknown")
        integrals = np.zeros(self.getOut())
        return integrals

