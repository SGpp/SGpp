// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractBasisCoefficientsStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>

#include <vector>
#include <memory>

namespace sgpp {
namespace combigrid {

class FunctionValuesCoefficientsStorage : public AbstractBasisCoefficientsStorage {
 public:
  FunctionValuesCoefficientsStorage();

  virtual ~FunctionValuesCoefficientsStorage();

  void computeCoefficients(MultiIndex const& level,
                           std::shared_ptr<AbstractCombigridStorage>& storage,
                           MultiIndex& multiBounds,
                           std::vector<bool>& orderingConfiguration) override;
};

} /* namespace combigrid */
} /* namespace sgpp */
