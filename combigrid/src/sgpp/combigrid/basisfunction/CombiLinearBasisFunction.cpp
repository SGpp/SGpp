// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include "CombiLinearBasisFunction.hpp"

using namespace std;

const combigrid::BasisFunctionBasis* combigrid::LinearBasisFunction::defaultBasis_ = new combigrid::LinearBasisFunction();
