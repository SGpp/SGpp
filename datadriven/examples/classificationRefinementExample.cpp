// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>

#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

/**
 * \page example_classificationRefinementExample_cpp Classification Example
 *
 * This example shows how classification specific refinement strategies
 * are used. To do classification, for each class a PDF is approximated with
 * LearnerSGDE and the class with the highest probability gets assigned for
 * new data points to be classified.
 * The ripley data sets is used, although the small number of training data
 * poitns in combination with only a basic setup does not yield good results
 * for any refinement strategy. This example is merely a tech-example.
 */

/**
 * Helper to create learner
 */
sgpp::datadriven::LearnerSGDE createSGDELearner(size_t dim, size_t level, double lambda);
/**
 * Helper to evaluate the classifiers
 */
std::vector<std::string> doClassification(std::vector<sgpp::base::Grid*> grids,
                                          std::vector<sgpp::base::DataVector*> alphas,
                                          sgpp::base::DataMatrix& testData,
                                          sgpp::base::DataVector& testLabel);

int main() {
  /**
   * Get the training/test data
   */
  std::string basePath = "../../datasets/ripley/ripleyGarcke";
  sgpp::datadriven::Dataset datasetTr =
    sgpp::datadriven::ARFFTools::readARFFFromFile(basePath + ".train.arff");
  sgpp::datadriven::Dataset datasetTs =
    sgpp::datadriven::ARFFTools::readARFFFromFile(basePath + ".test.arff");
  sgpp::base::DataMatrix dataTrain = datasetTr.getData();
  sgpp::base::DataVector targetTrain = datasetTr.getTargets();
  sgpp::base::DataMatrix dataTest = datasetTs.getData();
  sgpp::base::DataVector targetTest = datasetTs.getTargets();
  std::cout << "Read training data: " << dataTrain.getNrows() << std::endl;
  std::cout << "Read test data    : " << dataTest.getNrows() << std::endl;

  /**
   * Only uint class labels starting 0 and incrementing by 1 per class
   * are possible right no (to match grid-indices in vectors).
   * For use in DataVector class labels are cast to double.
   * Preprocess to have class label 0, 1, ...
   * -1 -> 0 and 1 -> 1
   */
  for (size_t i = 0; i < targetTrain.getSize(); i++) {
    if (targetTrain.get(i) < 0.0) {
      targetTrain.set(i, 0.0);
    } else {
      targetTrain.set(i, 1.0);
    }
  }
  for (size_t i = 0; i < targetTest.getSize(); i++) {
    if (targetTest.get(i) < 0.0) {
      targetTest.set(i, 0.0);
    } else {
      targetTest.set(i, 1.0);
    }
  }
  std::cout << "Preprocessing the data" << std::endl;

  /**
   * Split Training data according to class
   */
  sgpp::base::DataMatrix dataCl1(0.0, dataTrain.getNcols());
  sgpp::base::DataMatrix dataCl2(0.0, dataTrain.getNcols());
  sgpp::base::DataVector row(dataTrain.getNcols());
  for (size_t i = 0; i < dataTrain.getNrows(); i++) {
    dataTrain.getRow(i, row);
    if (targetTrain.get(i) < 1) {
      dataCl1.appendRow(row);
    } else {
      dataCl2.appendRow(row);
    }
  }
  std::cout << "Data points of class -1.0 (= 0): " << dataCl1.getNrows() << std::endl;
  std::cout << "Data points of class +1.0 (= 1): " << dataCl2.getNrows() << std::endl;

  /**
   * Approximate a probability density function for the class data using
   * LearnerSGDE, one for each class. Initialize the learners with the data
   */
  double lambda = 1e-5;
  sgpp::datadriven::LearnerSGDE learner1 = createSGDELearner(2, 3, lambda);
  sgpp::datadriven::LearnerSGDE learner2 = createSGDELearner(2, 3, lambda);
  learner1.initialize(dataCl1);
  learner2.initialize(dataCl2);

  /**
   * Bundle grids and surplus vector pointer needed for refinement and
   * evaluation
   */
  std::vector<sgpp::base::Grid*> grids;
  std::vector<sgpp::base::DataVector*> alphas;
  grids.push_back(learner1.getGrid());
  grids.push_back(learner2.getGrid());
  alphas.push_back(learner1.getSurpluses());
  alphas.push_back(learner2.getSurpluses());

  /**
   * Create refinement functors
   */
  size_t numRefinements = 3;
  bool levelPenalize = false;  // Multiplies penalzing term for fine levels
  bool preCompute = true;      // Precomputes and caches evals for zrcr & grid
  sgpp::datadriven::MultiGridRefinementFunctor* fun = nullptr;

  // Surplus refinement
  sgpp::datadriven::MultiSurplusRefinementFunctor funSrpl(grids, alphas, numRefinements,
                                                          levelPenalize);
  // Grid point-based refinement
  sgpp::datadriven::GridPointBasedRefinementFunctor funGrid(grids, alphas, numRefinements,
                                                            levelPenalize, preCompute);
  // Zero-crossing-based refinement
  sgpp::datadriven::ZeroCrossingRefinementFunctor funZrcr(grids, alphas, numRefinements,
                                                          levelPenalize, preCompute);
  /**
   * Data-based refinement. Needs a problem dependent coeffA. The values
   * were determined by testing (aim at ~10 % of the training data is
   * to be marked relevant. Cross-validation or similar can/should be employed
   * to determine this value.
   */
  std::vector<double> coeffA;
  coeffA.push_back(1.2);
  coeffA.push_back(1.2);
  sgpp::datadriven::DataBasedRefinementFunctor funData(grids, alphas, &dataTrain, &targetTrain,
                                                       numRefinements, levelPenalize, coeffA);

  /**
   * Choose the refinement functor to be used
   */
  // fun = &funSrpl;
  fun = &funGrid;
  // fun = &funZrcr;
  // fun = &funData;

  /**
   * Repeat alternating refinement and training for n steps and do
   * evaluation after each step
   * Uses the refinement strategy defined in fun
   * An initial evaluation with the regular grid is done at "step 0"
   */
  size_t numSteps = 5;
  std::vector<std::string> eval = doClassification(grids, alphas, dataTest, targetTest);
  std::cout << "Evaluation:" << std::endl << std::endl;
  std::cout << " Step  |  c=1    c=2  | total" << std::endl;
  std::cout << "------------------------------" << std::endl;
  std::cout << "   0   | " << eval.at(0) << " | " << eval.at(1) << std::endl;
  for (size_t i = 1; i < numSteps + 1; i++) {
    if (preCompute) {
      // precompute the evals. Needs to be done once per step, before
      // any refinement is done
      fun->preComputeEvaluations();
    }

    // Refine grid 0 for class -1.0 (= 0)
    fun->setGridIndex(0);
    grids.at(0)->getGenerator().refine(*fun);
    // Refine grid 1 for class 1.0 (= 1)
    fun->setGridIndex(1);
    grids.at(1)->getGenerator().refine(*fun);

    // Re-train the grids. This has to happend AFTER all grids are refined
    learner1.train(*grids.at(0), *alphas.at(0), dataCl1, lambda);
    learner2.train(*grids.at(1), *alphas.at(1), dataCl2, lambda);

    eval = doClassification(grids, alphas, dataTest, targetTest);
    std::cout << "   " << i << "   | " << eval.at(0) << " | " << eval.at(1) << std::endl;
  }
  std::cout << std::endl << "Done" << std::endl;
  return 0;
}

/**
 * Helper function
 * It configures and creates a SGDE learner with meaningful parameters
 */

sgpp::datadriven::LearnerSGDE createSGDELearner(size_t dim, size_t level, double lambda) {
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = dim;
  gridConfig.level_ = static_cast<int>(level);
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // configure adaptive refinement
  sgpp::base::AdaptivityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.noPoints_ = 10;

  // configure solver
  sgpp::solver::SLESolverConfiguration solverConfig;
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-10;

  // configure regularization
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Laplace;

  // configure learner
  sgpp::datadriven::CrossvalidationConfiguration crossvalidationConfig;
  crossvalidationConfig.enable_ = false;
  crossvalidationConfig.kfold_ = 3;
  crossvalidationConfig.lambda_ = 3.16228e-06;
  crossvalidationConfig.lambdaStart_ = lambda;
  crossvalidationConfig.lambdaEnd_ = lambda;
  crossvalidationConfig.lambdaSteps_ = 3;
  crossvalidationConfig.logScale_ = true;
  crossvalidationConfig.shuffle_ = true;
  crossvalidationConfig.seed_ = 1234567;
  crossvalidationConfig.silent_ = true;

  sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
                                        crossvalidationConfig);
  return learner;
}

/**
 * Helper function
 * it does the classification, gets the predictions and generates some error-output
 */

std::vector<std::string> doClassification(std::vector<sgpp::base::Grid*> grids,
                                          std::vector<sgpp::base::DataVector*> alphas,
                                          sgpp::base::DataMatrix& testData,
                                          sgpp::base::DataVector& testLabel) {
  double best_eval = 0.0;
  double eval = 0.0;
  sgpp::base::DataVector p(testData.getNcols());
  sgpp::base::DataVector indices(testData.getNrows());
  sgpp::base::DataVector evals(testData.getNrows());
  std::vector<std::unique_ptr<sgpp::base::OperationEval>> evalOps;
  for (size_t i = 0; i < grids.size(); i++) {
    std::unique_ptr<sgpp::base::OperationEval> e(
        sgpp::op_factory::createOperationEval(*grids.at(i)));
    evalOps.push_back(std::move(e));
  }

  // Get predictions and save to evals (confidence) and indices (class)
  for (size_t i = 0; i < testData.getNrows(); i++) {
    testData.getRow(i, p);
    indices.set(i, 0.0);
    best_eval = evalOps.at(0)->eval(*alphas.at(0), p);
    evals.set(i, best_eval);
    for (size_t j = 1; j < grids.size(); j++) {
      eval = evalOps.at(j)->eval(*alphas.at(j), p);
      if (eval > best_eval) {
        best_eval = eval;
        indices.set(i, static_cast<double>(j));
        evals.set(i, best_eval);
      }
    }
  }

  // Count the error for all classes
  sgpp::base::DataVector totalError(indices);
  std::vector<int> classCounts(grids.size(), 0);
  std::vector<int> classErrorCounts(grids.size(), 0);
  totalError.sub(testLabel);
  size_t totalCount = 0;
  for (size_t i = 0; i < testData.getNrows(); i++) {
    classCounts.at(static_cast<size_t>(floor(testLabel.get(i)))) += 1;
    if (fabs(totalError.get(i)) > 0.01) {
      totalCount++;
      classErrorCounts.at(static_cast<size_t>(floor(testLabel.get(i)))) += 1;
    }
  }

  // Format and return the classification percentages
  std::stringstream ss;
  for (size_t i = 0; i < grids.size(); i++) {
    double ce = (100.0 * (1.0 - (static_cast<double>(classErrorCounts.at(i)) / classCounts.at(i))));
    ss << std::fixed << std::setprecision(2) << ce;
    if (i < grids.size() - 1) {
      ss << "  ";
    }
  }
  std::stringstream ss2;
  ss2 << std::fixed << std::setprecision(3);
  ss2 << (100.0 *
          (1.0 - (static_cast<double>(totalCount) / static_cast<double>(testData.getNrows()))));
  std::vector<std::string> result;
  result.push_back(ss.str());
  result.push_back(ss2.str());
  return result;
}
