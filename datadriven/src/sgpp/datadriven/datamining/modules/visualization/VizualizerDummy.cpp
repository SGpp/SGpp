// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/visualization/VisualizerDummy.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>


using sgpp::datadriven::ModelFittingBase;

namespace sgpp {
namespace datadriven {

void VisualizerDummy::runVisualization(ModelFittingBase &model, DataSource &dataSource,
  size_t fold, size_t batch) {
  return;
}

}  // namespace datadriven
}  // namespace sgpp
