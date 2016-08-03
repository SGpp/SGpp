// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/MultiFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

class CombigridMultiOperationImpl;
// we use pimpl for not having to include all the template stuff
// in the header

class CombigridMultiOperation {
  std::shared_ptr<CombigridMultiOperationImpl> impl;  // unique_ptr causes SWIG errors

 public:
  CombigridMultiOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, MultiFunction func);

  CombigridMultiOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager,
      std::shared_ptr<AbstractCombigridStorage> storage);

  // TODO(holzmudd): add extra functions, for example for configuring the storage

  /**
   * Sets the parameters for upcoming computations and clears the data structures (removes old
   * computed data)
   */
  void setParameters(std::vector<base::DataVector> const &params);

  /**
   * Sets the parameters for upcoming computations and clears the data structures (removes old
   * computed data)
   */
  void setParameters(base::DataMatrix const &params);

  base::DataVector getResult();

  base::DataVector evaluate(
      size_t q,
      std::vector<base::DataVector> const &params);  // TODO(holzmudd): maybe change to DataMatrix?
  base::DataVector evaluateAdaptive(
      size_t maxNumPoints,
      std::vector<base::DataVector> const &params);  // TODO(holzmudd): maybe change to DataMatrix?

  std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> getDifferences();

  // TODO(holzmudd): add static constructor functions
  static std::shared_ptr<CombigridMultiOperation> createExpClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridMultiOperation> createLinearClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridMultiOperation> createExpLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridMultiOperation> createLinearLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);
  static std::shared_ptr<CombigridMultiOperation> createLinearLejaQuadrature(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);
  static std::shared_ptr<CombigridMultiOperation> createExpUniformLinearInterpolation(
      size_t numDimensions, MultiFunction func);
  // static std::shared_ptr<CombigridOperation> createLejaLinearPolynomialInterpolation(size_t
  // numDimensions); // TODO(holzmudd)
};
} /* namespace combigrid */
} /* namespace sgpp*/
