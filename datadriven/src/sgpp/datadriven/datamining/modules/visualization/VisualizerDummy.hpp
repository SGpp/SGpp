/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * VisualizerDummy.hpp
 *
 *  Created on: 16th Jun 2019
 *      Author: Vincent Bautista
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>

namespace sgpp {
namespace datadriven {

class VisualizerDummy:public Visualizer{

public:

 /**
  * Default constructor
  */
 VisualizerDummy()=default;


 ~VisualizerDummy()=default;


 void visualize(ModelFittingBase &model) override;

protected:

 void runTsne(ModelFittingBase &model) override;

};

}
}
