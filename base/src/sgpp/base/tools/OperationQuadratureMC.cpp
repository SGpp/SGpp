// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/OperationQuadratureMC.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstdlib>
#include <ctime>


namespace sgpp {
namespace base {


OperationQuadratureMC::OperationQuadratureMC(Grid& grid,
    int mcPaths) : grid(&grid), mcPaths(mcPaths) {
  // init seed for random number generator
  srand((unsigned)time(0));
}

double OperationQuadratureMC::doQuadrature(DataVector& alpha) {
  size_t dim = grid->getDimension();
  // create number of paths (uniformly drawn from [0,1]^d)
  DataMatrix dm(mcPaths, dim);

  for (size_t i = 0; i < mcPaths; i++) {
    for (size_t d = 0; d < dim; d++) {
      dm.set(i, d, static_cast<double>(rand()) / RAND_MAX);
    }
  }

  DataVector res = DataVector(mcPaths);
  sgpp::op_factory::createOperationMultipleEval(*grid, dm)->mult(alpha, res);
  return res.sum() / static_cast<double>(mcPaths);
}

double OperationQuadratureMC::doQuadratureFunc(FUNC func, void* clientdata) {
  size_t dim = grid->getDimension();
  double* p = new double[dim];

  // create number of paths (uniformly drawn from [0,1]^d)
  double res = 0;

  for (size_t i = 0; i < mcPaths; i++) {
    for (size_t d = 0; d < dim; d++) {
      p[d] = static_cast<double>(rand()) / RAND_MAX;
    }

    res += func(*reinterpret_cast<int*>(&dim), p, clientdata);
  }

  delete[] p;
  return res / static_cast<double>(mcPaths);
}

double OperationQuadratureMC::doQuadratureL2Error(FUNC func, void* clientdata,
    DataVector& alpha) {
  size_t dim = grid->getDimension();
  double x;
  double* p = new double[dim];

  sgpp::base::DataVector point(dim);
  std::unique_ptr<OperationEval> opEval = sgpp::op_factory::createOperationEval(*grid);
  // create number of paths (uniformly drawn from [0,1]^d)
  double res = 0;

  for (size_t i = 0; i < mcPaths; i++) {
    for (size_t d = 0; d < dim; d++) {
      x = static_cast<double>(rand()) / RAND_MAX;
      p[d] = x;
      point[d] = x;
    }

    res += pow(func(*reinterpret_cast<int*>(&dim), p,
                    clientdata) - opEval->eval(alpha, point), 2);
  }

  delete[] p;
  return sqrt(res / static_cast<double>(mcPaths));
}

}  // namespace base
}  // namespace sgpp
