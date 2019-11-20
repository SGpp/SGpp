// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

using sgpp::base::DataVector;

/**
 * We use metrics to quantify approximation quality of a trained model. It should be calculated in
 * the way that smaller is always better. This means, for example that classification accuracy
 * should be negated.
 */
class Metric {
 public:
  /**
   * Default constructor
   */
  Metric() = default;

  /**
   * Copy constructor
   * @param rhs const reference to the scorer object to copy from.
   */
  Metric(const Metric &rhs) = default;

  /**
   * Move constructor
   * @param rhs R-value reference to a scorer object to moved from.
   */
  Metric(Metric &&rhs) = default;

  /**
   * Copy assign operator
   * @param rhs const reference to the scorer object to copy from.
   * @return rerefernce to this with updated values.
   */
  Metric &operator=(const Metric &rhs) = default;

  /**
   * Move assign operator
   * @param rhs R-value reference to an a scorer object to move from.
   * @return rerefernce to this with updated values.
   */
  Metric &operator=(Metric &&rhs) = default;

  /**
   * virtual destructor.
   */
  virtual ~Metric() = default;

  /**
   * Polymorphic clone pattern
   * @return deep copy of this object. New object is owned by caller.
   */
  virtual Metric *clone() const = 0;

  /**
   * Quantify the difference between predicted values and actual values. Does not have an inner
   * state.
   *
   * @param predictedValues values calculated by the model for testing data
   * @param trueValues actual values as taken from the dataset.
   * @param model reference to the model
   * @param testDataset dataset with test data
   * @return Quantification of the difference. Smaller is better.
   */
  virtual double measure(const DataVector &predictedValues, const DataVector &trueValues,
                         const ModelFittingBase &model, Dataset &testDataset) const = 0;

};
} /* namespace datadriven */
} /* namespace sgpp */
