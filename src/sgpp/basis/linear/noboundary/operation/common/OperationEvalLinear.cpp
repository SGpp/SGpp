/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"

#include "basis/basis.hpp"

#include "basis/linear/noboundary/operation/common/OperationEvalLinear.hpp"

#include "data/DataVector.hpp"

namespace sg
{
namespace base
{

double OperationEvalLinear::eval(DataVector& alpha, std::vector<double>& point)
{
	linear_base<unsigned int, unsigned int> base;
	AlgorithmEvaluation<linear_base<unsigned int, unsigned int> > AlgoEval(storage);

	return AlgoEval(base, point, alpha);
}

}
}

