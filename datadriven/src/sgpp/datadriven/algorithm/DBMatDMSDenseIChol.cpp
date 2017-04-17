/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DBMatDMSDenseIChol.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: Michael Lettrich
 */

#include <sgpp/datadriven/algorithm/DBMatDMSDenseIChol.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/IChol.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

namespace sgpp {
namespace datadriven {
using sgpp::base::DataMatrix;
using sgpp::base::Grid;
using sgpp::base::OperationMatrix;

DBMatDMSDenseIChol::DBMatDMSDenseIChol(Grid& grid, double lambda, bool doCV)
    : DBMatDMSChol{}, proxyMatrix{} {
  // initialize proxy matrix if we do cv
  if (doCV) {
    auto size = grid.getStorage().getSize();
    proxyMatrix.resizeQuadratic(size);
    std::unique_ptr<OperationMatrix> op(
        sgpp::op_factory::createOperationLTwoDotExplicit(&proxyMatrix, grid));

    // set regularization parameter
    updateProxyMatrixLambda(lambda);
  }
}

void DBMatDMSDenseIChol::choleskyUpdateLambda(sgpp::base::DataMatrix& decompMatrix,
                                              double lambdaUpdate) const {
  updateProxyMatrixLambda(lambdaUpdate);
  IChol::decompose(proxyMatrix, decompMatrix, 2);
}

void DBMatDMSDenseIChol::updateProxyMatrixLambda(double lambdaUpdate) const {
  auto size = proxyMatrix.getNrows();
#pragma omp simd
  for (auto i = 0u; i < size; i++) {
    auto value = proxyMatrix.get(i, i) + lambdaUpdate;
    proxyMatrix.set(i, i, value);
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
