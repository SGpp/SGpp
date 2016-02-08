/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Janos Benk (benk@in.tum.de)



#include <sgpp/combigrid/basisfunction/CombiLinearBasisFunction.hpp>

using namespace std;

const combigrid::BasisFunctionBasis*
combigrid::LinearBasisFunction::defaultBasis_ = new
combigrid::LinearBasisFunction();

