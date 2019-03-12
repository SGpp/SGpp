// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/tools/MultipleClassPoint.hpp>

#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <cmath>

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

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
                                          sgpp::base::DataVector& testLabel, size_t classes);
/**
 * \page example_multipleClassExample_cpp Classification Example MultipleClassRefinement
 *
 * This example shows how the multiple class classification refinement strategy
 * is used. To do classification, for each class a PDF is approximated with
 * LearnerSGDE and the class with the highest probability gets assigned for
 * new data points to be classified.
 * This example is merely a tech-example.
 */

int main() {
  /**
   * All parameters are set in the beginning.
   * Allows to have an overview over set parameter.
   */
  // Parameter of data set
  std::string filepath = "../tests/data/";
  std::string filename = "multipleClassesTest.arff";
  // classes in ARFF are in [0,(classes-1)]
  size_t classes = 4;

  // Parameter for initial grid generation
  size_t dim = 2;
  size_t level = 4;
  double lambda = 1e-2;

  // Parameter for refinement
  size_t numSteps = 5;
  size_t numRefinements = 3;
  size_t partCombined = 0;
  double thresh = 0;

  // Only calculation after here, no additional parameters set

  sgpp::datadriven::Dataset dataset = sgpp::datadriven::ARFFTools::readARFFFromFile(filepath + filename);
  sgpp::base::DataMatrix dataTrain = dataset.getData();
  sgpp::base::DataVector targetTrain = dataset.getTargets();

  std::cout << "Read training data: " << dataTrain.getNrows() << std::endl;

  /**
   * Empty DataMartix are created to be filled with the data points from the data set
   * Using a vector, to be flexible for the amount of classes
   */
  std::vector<sgpp::base::DataMatrix> dataCl;
  std::vector<sgpp::datadriven::LearnerSGDE> learner;
  for (size_t i = 0; i < classes; i++) {
    dataCl.push_back(sgpp::base::DataMatrix(0.0, dataTrain.getNcols()));
  }

  /**
   * If classes are set to [0,classes-1] points are seperated into given classes.
   * Independent of the amount of classes needed
   *
   * Seperates the points into the different DataMatrix dependent on class
   */
  sgpp::base::DataVector row(dataTrain.getNcols());
  for (size_t i = 0; i < dataTrain.getNrows(); i++) {
    dataTrain.getRow(i, row);
    dataCl.at((size_t)targetTrain.get(i)).appendRow(row);
  }

  /**
   * Approximate a probability density function for the class data using
   * LearnerSGDE, one for each class. Initialize the learners with the data
   */
  for (size_t i = 0; i < classes; i++) {
    std::cout << "Data points of class " << std::setw(3) << std::right << i << ": ";
    std::cout << std::setw(14) << std::right << dataCl.at(i).getNrows() << "   | ";
    learner.push_back(createSGDELearner(dim, level, lambda));
    learner.back().initialize(dataCl.at(i));
  }

  /**
   * Bundle grids and surplus vector pointer needed for refinement and
   * evaluation
   */
  std::vector<sgpp::base::Grid*> grids;
  std::vector<sgpp::base::DataVector*> alphas;
  for (size_t i = 0; i < classes; i++) {
    grids.push_back(learner.at(i).getGrid());
    alphas.push_back(learner.at(i).getSurpluses());
  }

  std::cout << "---------------------------------------------" << std::endl;

  // Train the grids initially
  for (size_t i = 0; i < classes; i++) {
    learner.at(i).train(*grids.at(i), *alphas.at(i), dataCl.at(i), lambda);
  }
  std::vector<std::string> eval = doClassification(grids, alphas, dataTrain, targetTrain, classes);
  std::cout << "   0   | " << eval.at(0) << " | " << eval.at(1) << std::endl;

  // Create the MultipleClassRefinementFunctor with the set parameter
  sgpp::datadriven::MultipleClassRefinementFunctor mcrf(grids, alphas, numRefinements, partCombined,
                                                        thresh);
  sgpp::datadriven::MultipleClassRefinementFunctor* multifun = &mcrf;

  for (size_t x = 1; x < numSteps + 1; x++) {
    std::cout << "---------------------------------------------" << std::endl;

    // The refinement step is organizes by the Functor
    multifun->refine();

    // Re-train the grids. This has to happend AFTER all grids are refined
    for (size_t i = 0; i < classes; i++) {
      learner.at(i).train(*grids.at(i), *alphas.at(i), dataCl.at(i), lambda);
    }

    // Evaluate the result of the last refinement step
    eval = doClassification(grids, alphas, dataTrain, targetTrain, classes);
    std::cout << "   " << x << "   | " << eval.at(0) << " | " << eval.at(1) << std::endl;
  }
}

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
                                          sgpp::base::DataVector& testLabel, size_t classes) {
  double best_eval = -1000.0;
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
  std::vector<std::vector<int>> gridEval(classes + 1, std::vector<int>(classes + 1, 0));
  for (size_t i = 0; i < testData.getNrows(); i++) {
    testData.getRow(i, p);
    for (size_t j = 0; j < grids.size(); j++) {
      eval = evalOps.at(j)->eval(*alphas.at(j), p);
      if (eval > best_eval) {
        best_eval = eval;
        indices.set(i, static_cast<double>(j));
        evals.set(i, best_eval);
      }
    }
    best_eval = -1000.0;
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
