/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "algorithm/common/AlgorithmEvaluation.hpp"

#include "basis/basis.hpp"

#include "basis/linearstretched/noboundary/operation/common/OperationEvalLinearStretched.hpp"

#include "data/DataVector.hpp"

namespace sg
{
namespace base
{

double OperationEvalLinearStretched::eval(DataVector& alpha, std::vector<double>& point)
{
	linearstretched_base<unsigned int, unsigned int> base;
	AlgorithmEvaluation<linearstretched_base<unsigned int, unsigned int> > AlgoEval(storage);

	return AlgoEval(base, point, alpha);
}

}
}

