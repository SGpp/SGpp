// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>

#include <sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/data_exception.hpp>

namespace sgpp {

using sgpp::datadriven::NormalizationTransformationConfig;
using sgpp::datadriven::SampleProvider;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;

namespace datadriven {

class NormalizationTransformation : public DataTransformation {
 public:
  /**
   * Default constructor
   */

  NormalizationTransformation();

  /**
   * Initializes a transformation by determine minmax of each Dimension, if not given
   * @param dataset pointer to the dataset to be initialized
   * @param config configuration containing parameters for initalization
   */
  void initialize(Dataset *dataset, DataTransformationConfig config) override;

  /**
   * Wrapper for Normalization transformation. Can be called from an initialized
   * DataTransformation (with DataTransformationType::NORMALIZATION)
   *
   * @param dataset pointer to the dataset to be normalized
   * @return pointer to the transformed dataset
   */
  Dataset *doTransformation(Dataset *dataset) override;

  /**
   * Wrapper for Normalize inverse transformation (denormalized). Can be called from an initialized
   * DataTransformation (with DataTransformationType::NORMALIZATION)
   *
   * @param dataset pointer to the dataset to be denormalized
   * @return pointer to the backwards transformed dataset
   */
  Dataset *doInverseTransformation(Dataset *dataset) override;

 private:

  /**
   * Pointer to #sgpp::datadriven::NormalizationTransformationConfig
   */
  NormalizationTransformationConfig nmConfig;

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
