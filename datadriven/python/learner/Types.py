# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org



## Constants for types of grid borders
class BorderTypes(object):
    
    ## None border
    NONE = 100
    
    ## Trapezoid boundary
    TRAPEZOIDBOUNDARY = 200
    
    ## Complete boundary
    COMPLETEBOUNDARY = 300
    
    
## Constants for types of linear solvers    
class SolverTypes(object):
    
    ## Conjugate gradients (CG) solver
    CG = 100