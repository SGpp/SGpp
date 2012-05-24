/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonMLinear.hpp"

//#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
//#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"
//
//#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
//#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"
//
//#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
//#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"
//
//#include "pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp"

#include "base/algorithm/sweep.hpp"

#include <iostream>

namespace sg
{
namespace finance
{

OperationHestonMLinear::OperationHestonMLinear(sg::base::GridStorage* storage, double**** coef) : sg::pde::UpDownFourOpDims(storage, coef)
{
}

OperationHestonMLinear::~OperationHestonMLinear()
{
}

// Unidirectional
void OperationHestonMLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}

// Singles
void OperationHestonMLinear::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::downOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::upOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::downOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::upOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}

// Doubles
void OperationHestonMLinear::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::downOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::upOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::downOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::upOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::downOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::upOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonMLinear::downOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::upOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::downOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::upOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}

// Triples
void OperationHestonMLinear::downOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::upOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::downOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::upOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::downOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::upOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::downOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::upOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}

// Quadruples
void OperationHestonMLinear::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonMLinear::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}


}
}
