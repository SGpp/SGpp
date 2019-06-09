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
#include <sgpp/datadriven/datamining/modules/visualization/Parameters.hpp>


namespace sgpp {
namespace datadriven {

enum class Filetypes{csv,json};

enum class Libraries{plotly, mathplotlib, R, Matlab};

 struct VisualizationConfig{

  /**
   * The name of the algorithm to use in the visualization Module
   */
  std::string algorithm = "";

  /**
   * Parameters in which to run the algorithm
   */
  datadriven::Parameters;

  /**
   * The path to the file which will store the output of the visualization
   * module
   */
  std::string targetFile = "";


  /**
   * The filetype in which to store the output of the visualization module
   */
  std::string targetFileType = "";

  /**
   * The library to which the output will be forwarded
   */
  std::string targetLib = "";

};

}
}
