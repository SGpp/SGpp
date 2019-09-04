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

void Visualizer::createOutputDirectory(size_t fold, size_t batch) {
  std::cout << "Batch Number " << std::to_string(batch) << std::endl;
  if (config.getGeneralConfig().crossValidation) {
    currentDirectory = config.getGeneralConfig().
    targetDirectory+"/Fold_" + std::to_string(fold) + "/Batch_" + std::to_string(batch);
  } else {
    currentDirectory = config.getGeneralConfig().
    targetDirectory+"/Batch_" + std::to_string(batch);
  }

  std::cout << "Creating output directory " << config.getGeneralConfig().targetDirectory
     << std::endl;

  std::string mkdir("mkdir --parents ");

  mkdir.append(currentDirectory);

  int status = system(mkdir.data());

  if (status == 0) {
    std::cout << "Directory created succesfully" << std::endl;
  } else {
    std::cout << "Directory already exists" << std::endl;
  }
}

}  // namespace datadriven
}  // namespace sgpp
