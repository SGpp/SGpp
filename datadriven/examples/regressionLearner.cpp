#include <iostream>
#include <string.h>

#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include <sgpp/datadriven/tools/ARFFTools.hpp>

int main(int argc, char** argv) {

  //  int maxLevel = 9;
  int maxLevel = 5;

  //    std::string fileName = "debugging.arff";
  std::string fileName = "parabola_simple_3d.arff";
  //  std::string fileName = "friedman2_90000.arff";
  //  std::string fileName = "bigger.arff";

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset arffDataset = arffTools.readARFF(fileName);

  SGPP::base::DataMatrix& dataset = arffDataset.getTrainingData();
  SGPP::base::DataVector& values = arffDataset.getClasses();

  SGPP::base::RegularGridConfiguration gridConfig;
  SGPP::solver::SLESolverConfiguration SLESolverConfigRefine;
  SGPP::solver::SLESolverConfiguration SLESolverConfigFinal;
  SGPP::base::AdpativityConfiguration adaptConfig;

  // setup grid
  gridConfig.dim_ = arffDataset.getDimension(); //dim is inferred from the data
  gridConfig.level_ = maxLevel;
  gridConfig.type_ = SGPP::base::GridType::Linear;

  // Set Adaptivity
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  // Set solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 50;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = SGPP::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 50;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = SGPP::solver::SLESolverType::CG;

  std::string metaInformation = "refine: " + std::to_string((
                                  unsigned long long) adaptConfig.numRefinements_)
                                + " points: " + std::to_string((unsigned long long) adaptConfig.noPoints_) +
                                " iterations: "
                                + std::to_string((unsigned long long) SLESolverConfigRefine.maxIterations_);

  bool verbose = true;

  SGPP::datadriven::MetaLearner learner(gridConfig, SLESolverConfigRefine,
                                        SLESolverConfigFinal, adaptConfig,
                                        verbose);

  std::shared_ptr<SGPP::base::Grid> grid;
  std::shared_ptr<SGPP::base::DataVector> alpha;
  SGPP::float_t lambdaOpt;
  SGPP::datadriven::LearnerTiming timing;

  learner.optimizeLambdaLog(dataset, values, 3, 6, grid, alpha, lambdaOpt,
                            timing);

  SGPP::float_t mse = learner.calculateMSE(*grid, *alpha, dataset, values, true);
  std::cout << "final mse: " << mse << std::endl;

  return EXIT_SUCCESS;
}
