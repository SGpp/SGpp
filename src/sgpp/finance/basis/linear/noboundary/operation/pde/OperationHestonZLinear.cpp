/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonZLinear.hpp"

#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg
{
namespace finance
{

OperationHestonZLinear::OperationHestonZLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
{
}

OperationHestonZLinear::~OperationHestonZLinear()
{
}

void OperationHestonZLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
	sg::pde::PhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonZLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// phi * phi
	sg::pde::PhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonZLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiUpBBLinear func(this->storage);
	sg::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

void OperationHestonZLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{
	// x * dphi * phi
	XdPhiPhiDownBBLinear func(this->storage);
	sg::base::sweep<XdPhiPhiDownBBLinear> s(func, this->storage);

	s.sweep1D(alpha, result, dim);
}

}
}


///******************************************************************************
// * Copyright (C) 2009 Technische Universitaet Muenchen                         *
// * This file is part of the SG++ project. For conditions of distribution and   *
// * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
// ******************************************************************************/
//// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
//
//#include "finance/basis/linear/noboundary/operation/pde/OperationHestonZLinear.hpp"
//
//#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
//#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"
//
////#include "finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp"
////#include "finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp"
////
////#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
////#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"
////
////#include "finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinear.hpp"
////#include "finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinear.hpp"
//
////#include "finance/basis/linear/noboundary/operation/pde/OperationHestonXLinear.hpp"
////
//#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
//#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"
//
//#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
//#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"
//
//
//
//#include "base/algorithm/sweep.hpp"
//
//#include <iostream>
//
//namespace sg
//{
//namespace finance
//{
//
//OperationHestonZLinear::OperationHestonZLinear(sg::base::GridStorage* storage, sg::base::DataMatrix& coef) : sg::pde::UpDownTwoOpDims(storage, coef)
//{
//}
//
//OperationHestonZLinear::~OperationHestonZLinear()
//{
//}
//
//void OperationHestonZLinear::mult(sg::base::DataVector& alpha, sg::base::DataVector& result)
//{
//	result.setAll(0.0);
//
//#pragma omp parallel
//	{
//#pragma omp single nowait
//		{
//			for(size_t i = 0; i < this->numAlgoDims_; i++)
//			{
//				for(size_t j = 0; j < this->numAlgoDims_; j++)
//				{
//					// no symmetry in the operator
//#pragma omp task firstprivate(i, j) shared(alpha, result)
//					{
//						sg::base::DataVector beta(result.getSize());
//
//						if (this->coefs != NULL)
//						{
//							if (this->coefs->get(i,j) != 0.0)
//							{
//								this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j);
//
//#pragma omp critical
//								{
//									result.axpy(this->coefs->get(i,j),beta);
//								}
//							}
//						}
//						else
//						{
//							this->updown(alpha, beta, this->numAlgoDims_ - 1, i, j);
//
//#pragma omp critical
//							{
//								result.add(beta);
//							}
//						}
//					}
//				}
//			}
//
//#pragma omp taskwait
//		}
//	}
//}
//
//void OperationHestonZLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
//{
//	// phi * phi
//	sg::pde::PhiPhiUpBBLinear func(this->storage);
//	sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);
//
//	s.sweep1D(alpha, result, dim);
//}
//
//void OperationHestonZLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
//{
//	// phi * phi
//	sg::pde::PhiPhiDownBBLinear func(this->storage);
//	sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);
//
//	s.sweep1D(alpha, result, dim);
//}
//
//void OperationHestonZLinear::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
//{
//	// x * dphi * phi
//	XdPhiPhiUpBBLinear func(this->storage);
//	sg::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);
//
//	s.sweep1D(alpha, result, dim);
//}
//
//void OperationHestonZLinear::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
//{
//	// x * dphi * phi
//	sg::finance::XdPhiPhiDownBBLinear func(this->storage);
//	sg::base::sweep<sg::finance::XdPhiPhiDownBBLinear> s(func, this->storage);
//
//	s.sweep1D(alpha, result, dim);
//}
//
//void OperationHestonZLinear::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
//{
//	// phi * phi
//	sg::pde::PhiPhiUpBBLinear func(this->storage);
//	sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);
//
//	s.sweep1D(alpha, result, dim);
//}
//
//void OperationHestonZLinear::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
//{
//	// phi * phi
//	sg::pde::PhiPhiDownBBLinear func(this->storage);
//	sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);
//
//	s.sweep1D(alpha, result, dim);
//}
//
//void OperationHestonZLinear::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
//{
//}
//
//void OperationHestonZLinear::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
//{
//}
//
//}
//}
//
//
//
/////******************************************************************************
////* Copyright (C) 2009 Technische Universitaet Muenchen                         *
////* This file is part of the SG++ project. For conditions of distribution and   *
////* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
////******************************************************************************/
////// @author Sam Maurus (MA thesis)
////
////#include "finance/basis/linear/noboundary/operation/pde/OperationHestonXLinear.hpp"
////
////#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
////#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"
////
////#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
////#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"
////
////#include "base/algorithm/sweep.hpp"
////
////namespace sg
////{
////namespace finance
////{
////
////OperationHestonXLinear::OperationHestonXLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef)
////{
////}
////
////OperationHestonXLinear::~OperationHestonXLinear()
////{
////}
////
////void OperationHestonXLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
////{
////	// x * phi * phi
////	sg::finance::XPhiPhiUpBBLinear func(this->storage);
////	sg::base::sweep<sg::finance::XPhiPhiUpBBLinear> s(func, this->storage);
////
////	s.sweep1D(alpha, result, dim);
////}
////
////void OperationHestonXLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
////{
////	// x * phi * phi
////	sg::finance::XPhiPhiDownBBLinear func(this->storage);
////	sg::base::sweep<sg::finance::XPhiPhiDownBBLinear> s(func, this->storage);
////
////	s.sweep1D(alpha, result, dim);
////}
////
////void OperationHestonXLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
////{
////	// x * dphi * phi
////	sg::finance::XdPhiPhiUpBBLinear func(this->storage);
////	sg::base::sweep<sg::finance::XdPhiPhiUpBBLinear> s(func, this->storage);
////
////	s.sweep1D(alpha, result, dim);
////}
////
////void OperationHestonXLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
////{
////	// x * dphi * phi
////	sg::finance::XdPhiPhiDownBBLinear func(this->storage);
////	sg::base::sweep<sg::finance::XdPhiPhiDownBBLinear> s(func, this->storage);
////
////	s.sweep1D(alpha, result, dim);
////}
////
////}
////}
