from vectorFunctions import createJacobianErrorDataSet

model = 'demo'
path = '/home/rehmemk/git/SGpp/MR_Python/Vector/data/precalcGradients'
numMCPoints = 10000
dim = 2
out = 3
scalarModelParameter = 3  # dummy
createJacobianErrorDataSet(model, path, numMCPoints,
                           dim, out, scalarModelParameter)
