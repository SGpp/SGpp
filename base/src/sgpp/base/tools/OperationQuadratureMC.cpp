// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <cstdlib>
#include <ctime>

namespace sgpp {
namespace base {

OperationQuadratureMC::OperationQuadratureMC(Grid& grid, int mcPaths)
    : grid(&grid), mcPaths(mcPaths) {
  // init seed for random number generator
  srand(static_cast<unsigned>(time(nullptr)));
  // this->simple_rand.seed((unsigned)time(0));
}

double OperationQuadratureMC::doQuadrature(DataVector& alpha) {
  size_t dim = grid->getDimension();
  BoundingBox& boundingBox = grid->getBoundingBox();

  // create number of paths (uniformly drawn from [0,1]^d)
  DataMatrix dm(mcPaths, dim);

  for (size_t i = 0; i < mcPaths; i++) {
    for (size_t d = 0; d < dim; d++) {
      dm.set(i, d,
             boundingBox.transformPointToBoundingBox(d, static_cast<double>(rand()) / RAND_MAX));
      /*dm.set(i, d, boundingBox.transformPointToBoundingBox(
                       d, static_cast<double>(this->simple_rand()) / RAND_MAX));*/
    }
  }

  DataVector res(mcPaths);
  std::unique_ptr<OperationMultipleEval>(sgpp::op_factory::createOperationMultipleEval(*grid, dm))
      ->mult(alpha, res);

  // multiply with determinant of "unit cube -> BoundingBox" transformation
  double determinant = 1.0;

  for (size_t d = 0; d < dim; d++) {
    determinant *= boundingBox.getIntervalWidth(d);
  }

  return (res.sum() / static_cast<double>(mcPaths)) * determinant;
}

double OperationQuadratureMC::doQuadratureFunc(FUNC func, void* clientdata) {
  size_t dim = grid->getDimension();
  BoundingBox& boundingBox = grid->getBoundingBox();
  double* p = new double[dim];

  // create number of paths (uniformly drawn from [0,1]^d)
  double res = 0;

  for (size_t i = 0; i < mcPaths; i++) {
    for (size_t d = 0; d < dim; d++) {
      p[d] = boundingBox.transformPointToBoundingBox(d, static_cast<double>(rand()) / RAND_MAX);
      //p[d] = boundingBox.transformPointToBoundingBox(d, static_cast<double>(this->simple_rand()) / RAND_MAX);
    }

    res += func(static_cast<int>(dim), p, clientdata);
  }

  delete[] p;

  // multiply with determinant of "unit cube -> BoundingBox" transformation
  double determinant = 1.0;

  for (size_t d = 0; d < dim; d++) {
    determinant *= boundingBox.getIntervalWidth(d);
  }

  return res / static_cast<double>(mcPaths) * determinant;
}

double OperationQuadratureMC::doQuadratureL2Error(FUNC func, void* clientdata, DataVector& alpha) {
  size_t dim = grid->getDimension();
  BoundingBox& boundingBox = grid->getBoundingBox();
  double* p = new double[dim];

  sgpp::base::DataVector point(dim);
  std::unique_ptr<OperationEval> opEval(sgpp::op_factory::createOperationEval(*grid));
  // create number of paths (uniformly drawn from [0,1]^d)
  double res = 0;

  for (size_t i = 0; i < mcPaths; i++) {
    for (size_t d = 0; d < dim; d++) {
      p[d] = boundingBox.transformPointToBoundingBox(d, static_cast<double>(rand()) / RAND_MAX);
      //p[d] = boundingBox.transformPointToBoundingBox(d, static_cast<double>(this->simple_rand()) / RAND_MAX);
      point[d] = p[d];
    }

    res += pow(func(static_cast<int>(dim), p, clientdata) - opEval->eval(alpha, point), 2);
  }

  delete[] p;

  // multiply with determinant of "unit cube -> BoundingBox" transformation
  double determinant = 1.0;

  for (size_t d = 0; d < dim; d++) {
    determinant *= boundingBox.getIntervalWidth(d);
  }

  return sqrt(res / static_cast<double>(mcPaths) * determinant);
}

}  // namespace base
}  // namespace sgpp
