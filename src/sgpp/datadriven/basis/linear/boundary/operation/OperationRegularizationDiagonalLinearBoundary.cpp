/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "datadriven/basis/linear/boundary/operation/OperationRegularizationDiagonalLinearBoundary.hpp"
#include "base/grid/GridStorage.hpp"
#include <math.h>
//#include <iostream>

namespace sg
{
namespace datadriven
{

  OperationRegularizationDiagonalLinearBoundary::OperationRegularizationDiagonalLinearBoundary(base::GridStorage* storage, int mode, double k) 
      : OperationRegularizationDiagonal(storage, mode, k) {
      init();
    }

  void OperationRegularizationDiagonalLinearBoundary::initHkmix (double k) {
    int dim = storage->dim();
    base::GridIndex* gi;
    double res;
    for (size_t i=0; i<size; i++) {
      gi = storage->get(i);
      res = 1.0;
      for (int d=0; d<dim; d++) {
	res *= pow(2, (2*k-1) * gi->getLevel(d) -1);
      }
      diagonal[i] = res;
    }

  }
  
  void OperationRegularizationDiagonalLinearBoundary::initH0HkLaplace (double k) {
    int dim = storage->dim();
    base::GridIndex* gi;
    double res, resd;
    for (size_t i=0; i<size; i++) {
      gi = storage->get(i);
      res = 0.0;
      for (int d=0; d<dim; d++) {
	// Hk in dimension d
	resd = pow(2, (2*k-1) * gi->getLevel(d) -1);
	// "H0" in remaining dimensions
	for (int d2=0; d2<d; d2++) {
	  resd *= pow(2, -1 -gi->getLevel(d2));
	}
	for (int d2=d+1; d2<dim; d2++) {
	  resd *= pow(2, -1 -gi->getLevel(d2));
	}
	res += resd;
      }
      diagonal[i] = res;
    }

  }
  
}
}
