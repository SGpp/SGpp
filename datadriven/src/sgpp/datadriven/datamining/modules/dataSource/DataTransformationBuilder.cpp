// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformation.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/NormalizationTransformation.hpp>

namespace sgpp {
namespace datadriven {

DataTransformation *DataTransformationBuilder::buildTransformation(
    DataTransformationConfig config) {
  if (config.type == DataTransformationType::ROSENBLATT) {
    RosenblattTransformation *rosenblattTransformation = new RosenblattTransformation;
    return static_cast<DataTransformation *>(rosenblattTransformation);
  }
  if (config.type == DataTransformationType::NORMALIZATION) {
      NormalizationTransformation *normalizationTransformation = new NormalizationTransformation;
      return static_cast<DataTransformation *>(normalizationTransformation);
  }
  else {
    return nullptr;
  }
}
} /* namespace datadriven */
} /* namespace sgpp */
