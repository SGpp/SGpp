# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------
import anuga
import numpy as np
import shutil
from anuga.utilities import plot_utils as util
from anuga.config import netcdf_float
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from netCDF4 import Dataset
import numpy as np
from scipy.optimize import fmin


# par   input-wave parameters
# n     discretiazion size
def run():
    par = np.loadtxt('Okushiri/data/x.txt')
    n = np.loadtxt('Okushiri/data/gridsize.txt')
    # ------------------------------------------------------------------------------
    # Setup computational domain
    # ------------------------------------------------------------------------------
    # xleft = 0
    # xright = 5.448
    # ybottom = 0
    # ytop = 3.402
    # # rectangular cross mesh
    # points, vertices, boundary = anuga.rectangular_cross(int(n), int(n),
    #                                                      xright - xleft, ytop - ybottom,
    #                                                      (xleft, ybottom))

    # newpoints = points.copy()

    # # make refinement in x direction
    # x = np.multiply([0., 0.1, 0.2, 0.335, 0.925, 1.], max(points[:, 0]))
    # y = [0., 3., 4.25, 4.7, 5.3, max(points[:, 0])]
    # f1 = interp1d(x, y, kind='quadratic')
    # newpoints[:, 0] = f1(points[:, 0])

    # # make refinement in y direction
    # x = np.multiply([0., .125, .3, .7, .9, 1.], max(points[:, 1]))
    # y = [0., 1.25, 1.75, 2.15, 2.65, max(points[:, 1])]
    # f2 = interp1d(x, y, kind='quadratic')
    # newpoipar = np.loadtxt('Okushiri/data/x.txt')
    # n = np.loadtxt('Okushiri/data/gridsize.txt')
    # ------------------------------------------------------------------------------
    # Setup computational domain
    # ------------------------------------------------------------------------------
    xleft = 0
    xright = 5.448
    ybottom = 0
    ytop = 3.402

    # rectangular cross mesh
    points, vertices, boundary = anuga.rectangular_cross(int(n), int(n),
                                                         xright - xleft, ytop - ybottom,
                                                         (xleft, ybottom))

    newpoints = points.copy()

    # make refinement in x direction
    x = np.multiply([0., 0.1, 0.2, 0.335, 0.925, 1.], max(points[:, 0]))
    y = [0., 3., 4.25, 4.7, 5.3, max(points[:, 0])]
    f1 = interp1d(x, y, kind='quadratic')
    newpoints[:, 0] = f1(points[:, 0])

    # make refinement in y direction
    x = np.multiply([0., .125, .3, .7, .9, 1.], max(points[:, 1]))
    y = [0., 1.25, 1.75, 2.15, 2.65, max(points[:, 1])]
    f2 = interp1d(x, y, kind='quadratic')
    newpoints[:, 1] = f2(points[:, 1])

    c = abs(newpoints[:, 0] - 5.0) + .5 * abs(newpoints[:, 1] - 1.95)
    c = 0.125 * c

    points[:, 0] = c * points[:, 0] + (1 - c) * newpoints[:, 0]
    points[:, 1] = c * points[:, 1] + (1 - c) * newpoints[:, 1]

    # create domain
    domain = anuga.Domain(points, vertices, boundary)

    # don't store .sww file
    domain.set_quantities_to_be_stored(None)

    # ------------------------------------------------------------------------------
    # Initial Conditions
    # ------------------------------------------------------------------------------
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', 0.0)
    domain.set_quantity('elevation',
                        filename='/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/bathymetry.pts',
                        alpha=0.02)

    # ------------------------------------------------------------------------------
    # Set simulation parameters
    # ------------------------------------------------------------------------------
    domain.set_name('output_okushiri')  # Output name
    domain.set_minimum_storable_height(0.001)  # Don't store w < 0.001m
    domain.set_flow_algorithm('DE0')

    # ------------------------------------------------------------------------------
    # Modify input wave
    # ------------------------------------------------------------------------------
    # rescale input parameter
    try:
        dummy = len(par)
    except:
        par = [par]
    par = np.dot(2, par)

    # load wave data
    # shutil.copyfile('boundary_wave_header.txt', 'boundary_wave_input.txt')
    data = np.loadtxt(
        '/home/rehmemk/git/SGpp/MR_Python/Vector/Okushiri/boundary_wave_original.txt', skiprows=1)
    t = data[:, 0]
    y = data[:, 1]
    energy = np.trapz(y ** 2, t)

    # define bumps [create input wave based on parameters]
    def bump(c):
        theta = c[0]
        position = c[1]
        weight = c[2]
        ybump = weight * np.exp(-.5 * (t - position) ** 2 * theta ** -2)
        return ybump

    nbump = len(par)
    residual = y.copy()
    c = np.zeros((nbump, 3))
    for k in range(nbump):
        maxid = np.argmax(np.abs(residual))
        c0 = np.array([1.5, t[maxid], residual[maxid]])

        def cost(c):
            ybump = bump(c)
            cost = np.sqrt(np.mean((ybump - residual) ** 2))
            return cost

        c[k, :] = fmin(cost, c0, disp=False)
        residual -= bump(c[k, :])

    # deform wave
    ynew = residual.copy()
    for k in range(nbump):
        ynew += par[k] * bump(c[k, :])
    energynew = np.trapz(ynew ** 2, t)
    ynew = np.sqrt(energy / energynew) * ynew

    # write data
    data[:, 1] = ynew.copy()
    import scipy
    wave_function = scipy.interpolate.interp1d(
        t, ynew, kind='zero', fill_value='extrapolate')

    # MR: uncomment to plot input wave
    # points = np.linspace(-10, 30, 10000)
    # evals = np.zeros(len(points))
    # for i in range(len(evals)):
    #     evals[i] = wave_function(points[i])
    # plt.figure()
    # plt.plot(points, evals)
    # plt.plot(t, residual, 'r')
    # for k in range(nbump):
    #     plt.plot(t, par[k]*bump(c[k, :]))
    # plt.title('Okushiri Input Wave')
    # plt.show()

    # prepare time boundary
    # prepare_timeboundary('/home/rehmemk/git/SGpp/MR_Python/Extended/ANUGA/boundary_wave_input.tms')

    # ------------------------------------------------------------------------------
    # Setup boundary conditions
    # ------------------------------------------------------------------------------

    # Create boundary function from input wave [replaced by wave function]

    # Create and assign boundary objects
    Bts = anuga.Transmissive_momentum_set_stage_boundary(domain, wave_function)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'left': Bts, 'right': Br, 'top': Br, 'bottom': Br})

    # ------------------------------------------------------------------------------
    # Evolve system through time
    # ------------------------------------------------------------------------------

    # area for gulleys
    x1 = 4.85
    x2 = 5.25
    y1 = 2.05
    y2 = 1.85

    # index in gulley area
    x = domain.centroid_coordinates[:, 0]
    y = domain.centroid_coordinates[:, 1]
    v = np.sqrt((x - x1) ** 2 + (y - y1) ** 2) + \
        np.sqrt((x - x2) ** 2 + (y - y2) ** 2) < 0.5

    k = 0
    # original number of timesteps is 451
    numTimeSteps = int(np.loadtxt('Okushiri/data/numTimeSteps.txt'))
    meanstage = np.nan * np.ones((1, numTimeSteps))
    yieldstep = 0.05
    finaltime = (numTimeSteps - 1)*yieldstep
    for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
        # domain.write_time()

        # stage [=height of water]
        stage = domain.quantities['stage'].centroid_values[v]
        # averaging for smoothness
        meanstage[0, k] = np.mean(stage)
        # k is time
        k += 1

        # progress monitor
        np.savetxt('Okushiri/data/t.txt', [t])
    meanlayer = meanstage - meanstage[0, 0]
    # execnet can not serialze the ndarray meanlayer. Therefore save it to hard drive
    np.savetxt('Okushiri/data/y.txt', meanlayer)
    # return meanlayer


def prepare_timeboundary(filename, verbose=False):
    """Convert benchmark 2 time series to NetCDF tms file.
    This is a 'throw-away' code taylor made for files like
    'Benchmark_2_input.txt' from the LWRU2 benchmark
    """

    from anuga.file.netcdf import NetCDFFile

    if verbose:
        print 'Creating', filename

    # Read the ascii (.txt) version of this file
    fid = open(filename[:-4] + '.txt')

    # Skip first line
    line = fid.readline()

    # Read remaining lines
    lines = fid.readlines()
    fid.close()

    N = len(lines)
    T = np.zeros(N, np.float)  # Time
    Q = np.zeros(N, np.float)  # Values

    for i, line in enumerate(lines):
        fields = line.split()

        T[i] = float(fields[0])
        Q[i] = float(fields[1])

    # Create tms NetCDF file

    fid = NetCDFFile(filename, 'w')
    fid.institution = 'Geoscience Australia'
    fid.description = 'Input wave for Benchmark 2'
    fid.starttime = 0.0
    fid.createDimension('number_of_timesteps', len(T))
    fid.createVariable('time', netcdf_float, ('number_of_timesteps',))
    fid.variables['time'][:] = T

    fid.createVariable('stage', netcdf_float, ('number_of_timesteps',))
    fid.variables['stage'][:] = Q[:]

    fid.createVariable('xmomentum', netcdf_float, ('number_of_timesteps',))
    fid.variables['xmomentum'][:] = 0.0

    fid.createVariable('ymomentum', netcdf_float, ('number_of_timesteps',))
    fid.variables['ymomentum'][:] = 0.0

    fid.close()
