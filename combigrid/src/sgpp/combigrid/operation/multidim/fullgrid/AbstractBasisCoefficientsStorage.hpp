// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>

#include <vector>
#include <memory>
#include <map>

namespace sgpp {
namespace combigrid {

class AbstractBasisCoefficientsStorage {
 public:
  AbstractBasisCoefficientsStorage();

  virtual ~AbstractBasisCoefficientsStorage();
  virtual void computeCoefficients(MultiIndex const& level,
                                   std::shared_ptr<AbstractCombigridStorage>& storage,
                                   MultiIndex& multiBounds,
                                   std::vector<bool>& orderingConfiguration) = 0;
  std::shared_ptr<std::vector<double>> getCoefficients(
      MultiIndex const& level, std::shared_ptr<AbstractCombigridStorage>& storage,
      MultiIndex& multiBounds, std::vector<bool>& orderingConfiguration);

 protected:
  std::map<MultiIndex, std::shared_ptr<std::vector<double>>> coefficients;
};

} /* namespace combigrid */
} /* namespace sgpp */
