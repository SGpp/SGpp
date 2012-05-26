/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonKLinear.hpp"

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

OperationHestonKLinear::OperationHestonKLinear(sg::base::GridStorage* storage, double**** coef) : sg::pde::UpDownFourOpDims(storage, coef)
{
}

OperationHestonKLinear::~OperationHestonKLinear()
{
}

// Unidirectional
void OperationHestonKLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}

// Singles
void OperationHestonKLinear::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::downOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::upOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::downOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::upOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}

// Doubles
void OperationHestonKLinear::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::downOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::upOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::downOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::upOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::downOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::upOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim)
{

}
void OperationHestonKLinear::downOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}

// Triples
void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::downOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}

// Quadruples
void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}


}
}
