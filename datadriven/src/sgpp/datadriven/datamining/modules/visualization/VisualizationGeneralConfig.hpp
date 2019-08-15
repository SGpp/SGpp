/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * CrossvalidationConfiguration.hpp
 *
 *  Created on: Jun 8, 2019
 *      Author: Vincent Bautista
 */

#pragma once

#include <string>

namespace sgpp {
namespace datadriven {
enum class VisualizationFileType {CSV,json};

struct VisualizationGeneralConfig {

  /**
  * The name of the algorithm to use in the visualization Module
  */
  std::string algorithm = "";

  /**
  * The filetype in which to store the output of the visualization module
  */
  VisualizationFileType targetFileType =  VisualizationFileType::json;

  /**
  * Number of batches after which to execute the visualization Module
  * a 1 means every batch
  */
  size_t numBatches =  1;

  /**
  *  Flag that tells if cross validation is enabled. Depends on
  *  fitter[crossValidation][enabled]
  */
  bool crossValidation = false;

  /**
  * The path to the file which will store the output of the visualization
  * module
  */
  std::string targetDirectory = "";
};

}// namespace datadriven
}// namespace sgpp
