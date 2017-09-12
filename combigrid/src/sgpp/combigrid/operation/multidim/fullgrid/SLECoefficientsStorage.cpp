// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractBasisCoefficientsStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/SLECoefficientsStorage.hpp>

#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

SLECoefficientsStorage::SLECoefficientsStorage() : AbstractBasisCoefficientsStorage() {}

SLECoefficientsStorage::~SLECoefficientsStorage() {}

void SLECoefficientsStorage::computeCoefficients(MultiIndex const& level,
                                                 std::shared_ptr<AbstractCombigridStorage>& storage,
                                                 MultiIndex& multiBounds,
                                                 std::vector<bool>& orderingConfiguration) {
  MultiIndexIterator it(multiBounds);
  auto funcIter = storage->getGuidedIterator(level, it, orderingConfiguration);

  // compute interpolation matrix and collect function values
  base::DataVector functionValues;
  while (true) {
    // get function value and partial product and multiply them together with the last basis
    // coefficient, then add the resulting value to the total sum
    functionValues.push_back(funcIter->value());

    // increment iterator
    int h = funcIter->moveToNext();

    if (h < 0) {
      break;  // all indices have been traversed, stop iteration
    }
  }

  // ToDo (rehmemk) This is the SLE for Bspline interpolation only. Create a flag or something to
  // distinguish Bspline and potential other basis functions that might be used

  size_t numGridPoints = functionValues.size();
  base::DataMatrix A(numGridPoints, numGridPoints);
  A.setAll(0.0);
  for (size_t i = 0; i < numGridPoints; i++) {
    A.set(i, i, 1.0);
  }

  base::DataVector coefficients_sle(numGridPoints);

  optimization::Printer::getInstance().setVerbosity(-1);
  optimization::FullSLE sle(A);
  optimization::sle_solver::Armadillo solver;
  bool solved = solver.solve(sle, functionValues, coefficients_sle);

  if (!solved) {
    exit(-1);
  }

  coefficients[level] = std::make_shared<std::vector<double>>(coefficients_sle);
}

} /* namespace combigrid */
} /* namespace sgpp */
