/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataTransformation.cpp
 *
 *  Created on: 30.01.2018
 *      Author: Lars Wolfsteller
 */

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformation.hpp>

namespace sgpp {
namespace datadriven {

DataTransformation* DataTransformation::initialize(DataTransformationType dataTransformationType,
                                                  Dataset* dataset) {
  if (dataTransformationType == DataTransformationType::ROSENBLATT) {
    std::cout << "Initializing " <<  DataTransformationTypeParser::toString(dataTransformationType)
              << " transformation." << std::endl;
    return new RosenblattTransformation(dataset, 1000);
  } else {
    return new DataTransformation();
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
