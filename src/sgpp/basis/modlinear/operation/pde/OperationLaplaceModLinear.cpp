/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "basis/modlinear/operation/pde/OperationLaplaceModLinear.hpp"

#include "basis/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp"
#include "basis/modlinear/algorithm_sweep/dPhidPhiUpModLinear.hpp"
#include "basis/modlinear/algorithm_sweep/PhiPhiDownModLinear.hpp"
#include "basis/modlinear/algorithm_sweep/PhiPhiUpModLinear.hpp"

#include "algorithm/common/sweep.hpp"
using namespace sg::base;

namespace sg
{
namespace pde
{

OperationLaplaceModLinear::OperationLaplaceModLinear(GridStorage* storage) : UpDownOneOpDim(storage)
{
}

OperationLaplaceModLinear::~OperationLaplaceModLinear()
{
}

void OperationLaplaceModLinear::up(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	PhiPhiUpModLinear func(this->storage);
	sweep<PhiPhiUpModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	PhiPhiDownModLinear func(this->storage);
	sweep<PhiPhiDownModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	dPhidPhiDownModLinear func(this->storage);
	sweep<dPhidPhiDownModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	dPhidPhiUpModLinear func(this->storage);
	sweep<dPhidPhiUpModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

}
}

