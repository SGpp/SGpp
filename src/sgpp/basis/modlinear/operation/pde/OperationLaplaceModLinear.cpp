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

namespace sg
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
	detail::PhiPhiUpModLinear func(this->storage);
	sweep<detail::PhiPhiUpModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::down(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	detail::PhiPhiDownModLinear func(this->storage);
	sweep<detail::PhiPhiDownModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::downOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	detail::dPhidPhiDownModLinear func(this->storage);
	sweep<detail::dPhidPhiDownModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

void OperationLaplaceModLinear::upOpDim(DataVector& alpha, DataVector& result, size_t dim)
{
	result.setAll(0.0);
	detail::dPhidPhiUpModLinear func(this->storage);
	sweep<detail::dPhidPhiUpModLinear> s(func, this->storage);
	s.sweep1D(alpha, result, dim);
}

}

