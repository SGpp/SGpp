/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "common/operation/OperationQuadratureMC.hpp"
#include "basis/operations_factory.hpp"
#include "data/DataMatrix.hpp"

namespace sg
{
namespace base
{
  OperationQuadratureMC::OperationQuadratureMC(Grid &grid, int mcPaths) : grid(&grid), mcPaths(mcPaths) {
    // init seed for random number generator
    srand((unsigned)time(0)); 
  }


double OperationQuadratureMC::doQuadrature(DataVector& alpha)
{
  int dim = grid->getStorage()->dim();
  // create number of paths (uniformly drawn from [0,1]^d)
  DataMatrix dm(mcPaths, dim);
  for (int i=0; i<mcPaths; i++) {
    for (int d=0; d<dim; d++) {
      dm.set(i, d, static_cast<double>(rand())/RAND_MAX);
    }
  }
  OperationMultipleEval* opEval = sg::GridOperationFactory::createOperationMultipleEval(*grid, &dm);
  DataVector res = DataVector(mcPaths);
  opEval->mult(alpha, res);
  return res.sum()/static_cast<double>(mcPaths);
}

}
}
