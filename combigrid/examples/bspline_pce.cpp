// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModelFactory.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/combigrid/utils/AnalyticModels.hpp>

#include <cmath>

#include <iostream>
#include <vector>

double u(sgpp::base::DataVector x) {
  return 4. * x[0] * (1. - x[0]);
  //  return std::log(std::exp(-x[0]) * std::cos(4. * x[0] * (1. - x[0])));
}

int main() {
  try {
    sgpp::combigrid::Parabola_uniform model;
    sgpp::combigrid::MultiFunction func(model.eval);

    size_t numDims = 2;
    size_t q = 1;
    size_t degree = 3;

    // interpolate on adaptively refined grid
    auto op = sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineInterpolation(
        numDims, func, degree);
    op->getLevelManager()->addRegularLevels(q);

    // compute mean
    auto quad_op = sgpp::combigrid::CombigridOperation::createExpUniformBoundaryBsplineQuadrature(
        numDims, func, degree);
    quad_op->getLevelManager()->addLevelsFromStructure(op->getLevelManager()->getLevelStructure());
    double mean = quad_op->getResult();

    auto tensor_op =
        sgpp::combigrid::CombigridTensorOperation::createExpUniformBoundaryBSplineInterpolation(
            numDims, func, degree);
    tensor_op->getLevelManager()->addLevelsFromStructure(
        op->getLevelManager()->getLevelStructure());

    sgpp::combigrid::FloatTensorVector tensor_result = tensor_op->getResult();

    size_t numTerms = 0;
    double variance = 0.0;
    std::cout << " ----------------------------------" << std::endl;
    for (auto it = tensor_result.getValues()->getStoredDataIterator(); it->isValid();
        it->moveToNext()) {
      double ocoeff = it->value().value();
      variance += ocoeff * ocoeff;
      numTerms += 1;
      auto ix = it->getMultiIndex();
      std::cout << numTerms << ": (" << ix[0] << ", " << ix[1] << ") -> " << ocoeff << std::endl;
    }

    variance -= mean * mean;

    std::cout << "---------------------------------------------------------" << std::endl;
    std::cout << "#gp    = " << tensor_op->getLevelManager()->numGridPoints() << std::endl;
    std::cout << "#terms = " << numTerms << std::endl;
    std::cout << "E(u)   = " << model.mean(numDims) << " ~ " << mean << std::endl;
    std::cout << "Var(u) = " << model.variance(numDims) << " ~ " << variance << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;
  }
  catch (sgpp::base::generation_exception& exc)  {
    std::cout << "Exception: " << exc.what() << std::endl;
    std::cout << "Skipping example..." << std::endl;
  }
}
