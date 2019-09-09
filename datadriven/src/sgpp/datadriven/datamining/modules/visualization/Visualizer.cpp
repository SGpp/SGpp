/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Visualizer.hpp
 *
 *  Created on: 16th Jun 2019
 *      Author: Vincent Bautista
 */

#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <string>

namespace sgpp {
namespace datadriven {

Visualizer::Visualizer() {
}

const VisualizerConfiguration &Visualizer::getVisualizerConfiguration() const {
  return config;
}

void Visualizer::createFolder(std::string folder_path) {
  #ifdef _WIN32
  std::string mkdir("mkdir -p ");
  #elif __linux__
  std::string mkdir("mkdir --parents ");
  #endif

  mkdir.append(folder_path);

  int status = system(mkdir.data());

  if (status == 0) {
    std::cout << "Folder creation command executed succesfully" << std::endl;
  } else {
    std::cout << "Error executing folder creation command" << std::endl;
  }

}

}  // namespace datadriven
}  // namespace sgpp
