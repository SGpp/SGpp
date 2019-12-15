// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {

Scorer::Scorer(Metric* metric) : metric{std::unique_ptr<Metric>{metric}} {}

double Scorer::test(ModelFittingBase& model, Dataset& testDataset, bool lowerIsBetter) {
#ifdef USE_SCALAPACK
  if (model.getFitterConfiguration().getParallelConfig().scalapackEnabled_) {
    return testDistributed(model, testDataset, lowerIsBetter);
  }
#endif
  DataVector predictedValues(testDataset.getNumberInstances());
  model.evaluate(testDataset.getData(), predictedValues);
  // set score
  if (lowerIsBetter) {
    return metric->measureLowerIsBetter(predictedValues, testDataset.getTargets(), model,
                                        testDataset);
  }
  return metric->measure(predictedValues, testDataset.getTargets(), model, testDataset);
}

double Scorer::testDistributed(ModelFittingBase& model, Dataset& testDataset, bool lowerIsBetter) {
#ifdef USE_SCALAPACK
  DataVector predictedValues(testDataset.getNumberInstances());
  model.evaluate(testDataset.getData(), predictedValues);

  // only calculate the score on the master and send the result to the other processes, this means
  // the evaluation results don't have to broadcast.
  double score = 0.0;
  auto processGrid = model.getProcessGrid();
  if (processGrid->getCurrentRow() == 0 && processGrid->getCurrentColumn() == 0) {
    if (lowerIsBetter) {
      score = metric->measureLowerIsBetter(predictedValues, testDataset.getTargets(), model,
                                           testDataset);
    } else {
      score = metric->measure(predictedValues, testDataset.getTargets(), model, testDataset);
    }
    Cdgebs2d(processGrid->getContextHandle(), "All", "T", 1, 1, &score, 1);
  } else if (processGrid->isProcessInGrid()) {
    Cdgebr2d(processGrid->getContextHandle(), "All", "T", 1, 1, &score, 1, 0, 0);
  } else {
    std::cout << "Warning! Process not in the grid tried to call testDistributed, invalid result "
                 "will be returned!"
              << std::endl;
  }

  return score;
#else
  throw sgpp::base::application_exception("Build without USE_SCALAPACK");
#endif /* USE_SCALAPACK */
}

} /* namespace datadriven */
}  // namespace sgpp
