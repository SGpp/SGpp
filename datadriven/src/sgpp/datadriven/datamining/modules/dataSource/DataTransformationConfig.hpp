// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformationConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/NormalizationTransformationConfig.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Supported transformation types for sgpp::datadriven::DataTransformation
 */
enum class DataTransformationType { NONE, ROSENBLATT, NORMALIZATION };

struct DataTransformationConfig {
  /*
   * Type of data transformation
   */
  DataTransformationType type = DataTransformationType::NONE;

  RosenblattTransformationConfig rosenblattConfig;

  NormalizationTransformationConfig normalizationConfig;
};
} /* namespace datadriven */
} /* namespace sgpp */
