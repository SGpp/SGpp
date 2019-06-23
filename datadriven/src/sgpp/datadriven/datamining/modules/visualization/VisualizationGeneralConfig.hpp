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

#include "VisualizationParameters.hpp"


namespace sgpp {
namespace datadriven {

enum class VisualizationFileType{CSV,json};


 struct VisualizationGeneralConfig{

  /**
   * The name of the algorithm to use in the visualization Module
   */
  std::string algorithm = "";

  /**
   * The path to the file which will store the output of the visualization
   * module
   */
  std::string targetFile ="";


  /**
   * The filetype in which to store the output of the visualization module
   */
  VisualizationFileType targetFileType =  VisualizationFileType::json;


};

}
}
