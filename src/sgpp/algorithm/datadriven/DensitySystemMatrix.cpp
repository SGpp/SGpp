/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Fabian Franzelin <franzeli@in.tum.de>,
//         Benjamin Peherstorfer <pehersto@in.tum.de>

#include "algorithm/datadriven/DensitySystemMatrix.hpp"
#include "basis/operations_factory.hpp"
#include "exception/operation_exception.hpp"
#include "basis/linear/noboundary/operation/pde/OperationLTwoDotProductLinear.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

namespace sg
{
namespace datadriven
{

  DensitySystemMatrix::DensitySystemMatrix(sg::base::Grid& grid, sg::base::DataMatrix& trainData, sg::base::OperationMatrix& C, double lambda)
{
  this->data = &trainData;
  this->lambda = lambda;

  this->A = sg::op_factory::createOperationLTwoDotProduct(grid);
  this->B = sg::op_factory::createOperationMultipleEval(grid, this->data);
  this->C = &C;
}

void DensitySystemMatrix::mult(sg::base::DataVector &alpha, sg::base::DataVector &result)
{
  result.setAll(0.0);

  // A * alpha
  this->A->mult(alpha, result);
  
  // C * alpha
  base::DataVector tmp(result.getSize());
  this->C->mult(alpha, tmp);

  // A * alpha + lambda * C * alpha
  result.axpy(this->lambda, tmp);
}


// Matrix-Multiplikation verwenden
void DensitySystemMatrix::generateb(sg::base::DataVector& rhs)
{
  sg::base::DataVector y(this->data->getNrows());
  y.setAll(1.0);
  // Bt * 1
  this->B->multTranspose(y, rhs);
  // 1 / 2M * Bt * 1
  rhs.mult(1./ (2. * this->data->getNrows()));
}

DensitySystemMatrix::~DensitySystemMatrix()
{
  delete this->A;
  delete this->B;
}

}
}
