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
#include <string>
#include<iostream>
#ifdef _WIN32
#include <direct.h>
#elif defined __linux__
#include <sys/stat.h>
#include <sys/types.h>
#endif
namespace sgpp {
namespace datadriven {

Visualizer::Visualizer() {
}

const VisualizerConfiguration &Visualizer::getVisualizerConfiguration() const {
  return config;
}

void Visualizer::createFolder(std::string folder_path) {
  #ifdef _WIN32
  _mkdir(folder_path.c_str());
  #elif __linux__
  mkdir(folder_path.c_str(), S_IRWXU | S_IRWXU | S_IROTH);
  #endif
}

}  // namespace datadriven
}  // namespace sgpp
