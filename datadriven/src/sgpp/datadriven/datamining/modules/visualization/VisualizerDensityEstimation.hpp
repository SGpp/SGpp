/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * VisualizerDensityEstimation.hpp
 *
 *  Created on: 16th Jun 2019
 *      Author: Vincent Bautista
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizerConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>
#include <sgpp/datadriven/tools/CSVTools.hpp>
#include <vector>
namespace sgpp {
namespace datadriven {

class VisualizerDensityEstimation:public Visualizer{

public:

 /**
  * Default constructor
  */
 VisualizerDensityEstimation()=default;

 /**
  * Constructor given a configuration
  * @param config. The VisualizerConfiguration object which contains
  * the configuration to run the visualization module
  */
 VisualizerDensityEstimation(VisualizerConfiguration config);

 ~VisualizerDensityEstimation()=default;


 void visualize(ModelFittingBase &model) override;

protected:

 void run_tsne() override;

 void getHeatmap(ModelFittingBase &model);

 void getLinearCuts(ModelFittingBase &model);

private:
 void translateColumns(DataMatrix &matrix, size_t maxColumns);
 void translateColumns(DataMatrix &matrix, std::vector<size_t> indexes);
 void updateIndexes(std::vector<size_t> &columnIndexes);


};

}
}
