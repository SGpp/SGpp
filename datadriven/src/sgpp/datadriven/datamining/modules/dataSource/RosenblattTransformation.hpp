// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

namespace sgpp {

using sgpp::base::Grid;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::datadriven::LearnerSGDE;

namespace datadriven {

class RosenblattTransformation : public DataTransformation {
 public:
  /**
   * Default constructor
   */

  RosenblattTransformation();

  /**
   * Initializes a transformation by approximating probability density function (PDF),
   * calculates grid and alpha for numSamples samples of a dataset
   * @param dataset pointer to the dataset to be initialized
   * @param config configuration containing parameters for initalization
   */
  void initialize(Dataset *dataset, DataTransformationConfig config) override;

  /**
   * Wrapper for Rosenblatt transformation. Can be called from an initialized
   * DataTransformation (with DataTransformationType::ROSENBLATT)
   *
   * @param dataset pointer to the dataset to be Rosenblatt transformed
   * @return pointer to the transformed dataset
   */
  Dataset *doTransformation(Dataset *dataset) override;

  /**
   * Wrapper for Rosenblatt inverse transformation. Can be called from an initialized
   * DataTransformation (with DataTransformationType::ROSENBLATT)
   *
   * @param dataset pointer to the dataset to be Rosenblatt inverse transformed
   * @return pointer to the backwards transformed dataset
   */
  Dataset *doInverseTransformation(Dataset *dataset) override;

  /**
   * Helper function
   * It configures and creates a SGDE learner with meaningful parameters
   */
  LearnerSGDE createSGDELearner(size_t dim, RosenblattTransformationConfig config);

 private:
  /**
     * the sparse grid that approximates the data.
     */
  std::shared_ptr<base::Grid> grid;

  /**
   * hierarchical surpluses of the #grid.
   */
  std::shared_ptr<base::DataVector> alpha;

  /**
   * Pointer to #sgpp::datadriven::Dataset
   */
  Dataset *datasetTransformed;

  /**
   * Pointer to #sgpp::datadriven::Dataset
   */
  Dataset *datasetInvTransformed;
};
} /* namespace datadriven */
} /* namespace sgpp */
