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

#pragma once

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>

namespace sgpp {

namespace datadriven {

class Visualizer{

 public:

 /*
  * Default Constructor
  */
  Visualizer()=default;

 /*
  * Virtual destructor
  */
  virtual ~Visualizer()=default;

 /**
  * Method to execute the visualization step
  */
  virtual void visualize(ModelFittingBase &model)=0;

 /**
   * Get the configuration of the visualizer object.
   * @return configuration of the visualizer object
   */
   const VisualizerConfiguration &getVisualizerConfiguration() const;


 protected:

  virtual void run_tsne(ModelFittingBase &model)=0;

  void createOutputDirectory();


  /**
    * Configuration object for the fitter.
    */
  VisualizerConfiguration config;


};
}



}
