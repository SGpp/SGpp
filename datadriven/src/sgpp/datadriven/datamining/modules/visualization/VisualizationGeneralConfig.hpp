// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#pragma once

#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {
enum class VisualizationFileType {CSV, json};

struct VisualizationGeneralConfig {
  /**
   * Variable to determine if its executes visualization Module or not
   */
  bool execute_ = false;

  /**
  * The list of algorithms to use in the visualization Module
  */
  std::vector<std::string> algorithm_ = std::vector<std::string>();

  /**
  * The filetype in which to store the output of the visualization module
  */
  VisualizationFileType targetFileType_ =  VisualizationFileType::json;

  /**
  * Number of batches after which to execute the visualization Module
  * a 1 means every batch
  */
  size_t numBatches_ =  1;

  /**
  *  Flag that tells if cross validation is enabled. Depends on
  *  fitter[crossValidation][enabled]
  */
  bool crossValidation_ = false;

  /**
  * The path to the file which will store the output of the visualization
  * module
  */
  std::string targetDirectory_ = "";
};
}  // namespace datadriven
}  // namespace sgpp
