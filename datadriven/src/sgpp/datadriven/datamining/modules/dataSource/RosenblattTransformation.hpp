/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * RosenblattTransformation.hpp
 *
 *  Created on: 22.01.2018
 *      Author: Lars Wolfsteller
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

namespace sgpp {

using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;

namespace datadriven {

class RosenblattTransformation : public DataTransformation {
 public:
  /**
   * Default constructor
   */
  RosenblattTransformation(Dataset* dataset, size_t numSamples);

  /**
   * Wrapper for Rosenblatt transformation. Can be called from an initialized
   * DataTransformation (with DataTransformationType::ROSENBLATT)
   *
   * @param dataset pointer to the dataset to be Rosenblatt transformed
   * @return pointer to the transformed dataset
   */
  Dataset* doTransformation(Dataset* dataset) override;

  /**
   * Wrapper for Rosenblatt inverse transformation. Can be called from an initialized
   * DataTransformation (with DataTransformationType::ROSENBLATT)
   *
   * @param dataset pointer to the dataset to be Rosenblatt inverse transformed
   * @return pointer to the backwards transformed dataset
   */
  Dataset* doInverseTransformation(Dataset* dataset) override;

  /**
   * Helper function
   * It configures and creates a SGDE learner with meaningful parameters
   */
  sgpp::datadriven::LearnerSGDE createSGDELearner(size_t dim, size_t level,
                                                    double lambda);

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
  Dataset* datasetTransformed;

  /**
   * Pointer to #sgpp::datadriven::Dataset
   */
  Dataset* datasetInvTransformed;

  /**
   * Number of samples for calculation of pdf / alpha
   */
  size_t numSamples;
};
} /* namespace datadriven */
} /* namespace sgpp */
