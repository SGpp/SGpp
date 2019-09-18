from vectorFunctions import createJacobianErrorDataSet

model = 'dc_motor_I'
path = '/home/rehmemk/git/SGpp/MR_Python/Vector/data/precalcGradients'
numMCPoints = 1000
dim = 6
out = 101
scalarModelParameter = 3  # dummy
createJacobianErrorDataSet(model, path, numMCPoints,
                           dim, out, scalarModelParameter)
