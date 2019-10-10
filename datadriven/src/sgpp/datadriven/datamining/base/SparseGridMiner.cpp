// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/algorithm/RefinementMonitorFactory.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <iostream>
#include <string>

namespace sgpp {
namespace datadriven {
SparseGridMiner::SparseGridMiner(ModelFittingBase* fitter, Scorer* scorer, Visualizer *visualizer)
    : fitter(fitter), scorer(scorer), visualizer(visualizer) {}

double SparseGridMiner::test(Dataset& testDataset) { return scorer->test(*fitter, testDataset); }

ModelFittingBase* SparseGridMiner::getModel() { return &(*fitter); }

void SparseGridMiner::setModel(ModelFittingBase* model) { fitter.reset(model); }

Visualizer* SparseGridMiner::getVisualizer() { return &(*visualizer); }

void SparseGridMiner::print(std::ostringstream& messageStream) { print(messageStream.str()); }

void SparseGridMiner::print(const char* message) { print(std::string(message)); }

void SparseGridMiner::print(const std::string& message) {
#ifdef USE_SCALAPACK
  if (BlacsProcessGrid::getCurrentProcess() == 0) {
#endif
    std::cout << message << std::endl;
#ifdef USE_SCALAPACK
  }
#endif
}

} /* namespace datadriven */
}  // namespace sgpp
