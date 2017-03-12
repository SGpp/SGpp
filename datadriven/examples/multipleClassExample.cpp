// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>


#include <cmath>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sgpp/datadriven/application/MultipleClassPoint.hpp>


/**
 * Helper to create learner
 */
sgpp::datadriven::LearnerSGDE createSGDELearner(size_t dim, size_t level,
                                                double lambda);
 
/**
 * Helper to evaluate the classifiers
 */
std::vector<std::string> doClassification(std::vector<sgpp::base::Grid*> grids,
                                          std::vector<sgpp::base::DataVector*> alphas,
                                          sgpp::base::DataMatrix& testData,
                                          sgpp::base::DataVector& testLabel);
/**
 * TODO
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


int main() {

    std::string filename2 = "../../datadriven/tests/data/mulitpleClassesTest.arff";
    sgpp::datadriven::Dataset dataset =
            sgpp::datadriven::ARFFTools::readARFF(filename2);
    sgpp::base::DataMatrix dataTrain = dataset.getData();
    sgpp::base::DataVector targetTrain = dataset.getTargets();
    std::cout << "Read training data: " << dataTrain.getNrows() << std::endl;
    int classes = 4;

    std::vector<sgpp::base::DataMatrix> dataCl;
    std::vector<sgpp::datadriven::LearnerSGDE> learner;
    for ( int i = 0 ; i < classes ; i++ ) {
        dataCl.push_back(sgpp::base::DataMatrix(0.0, dataTrain.getNcols()));
    }

    sgpp::base::DataVector row(dataTrain.getNcols());
    for ( size_t i = 0 ; i < dataTrain.getNrows() ; i++ ) {
        dataTrain.getRow(i, row);
        dataCl.at((int)targetTrain.get(i)).appendRow(row);
    }
    for ( int i = 0 ; i < classes ; i++ ) {
        std::cout << "Data points of class " << i << ": ";
        std::cout << dataCl.at(i).getNrows() << std::endl;
    }

    /**
    * Approximate a probability density function for the class data using
    * LearnerSGDE, one for each class. Initialize the learners with the data
    */
    size_t dim = 2;
    int level = 3;
    double lambda = 1e-1;//größer (start 5)
    for ( int i = 0 ; i < classes ; i++ ) {
        learner.push_back(createSGDELearner(dim, level, lambda));
        learner.back().initialize(dataCl.at(i));
    }

    /**
     * Bundle grids and surplus vector pointer needed for refinement and
     * evaluation
     */
    std::vector<sgpp::base::Grid*> grids2;
    std::vector<sgpp::base::DataVector*> alphas2;
    for ( int i = 0 ; i < classes ; i++ ) {
        grids2.push_back(learner.at(i).getGrid().get());
        alphas2.push_back(learner.at(i).getSurpluses().get());
    }
    for ( int i = 0 ; i < classes ; i++ ) {
        learner.at(i).train(*grids2.at(i), *alphas2.at(i), dataCl.at(i), lambda);
    }

    size_t numSteps = 1;
    std::vector<std::string> eval = doClassification(grids2, alphas2, dataTrain, targetTrain);

    std::cout << std::endl << "Start" << std::endl;
    size_t numRefinements = 3;
    bool levelPenalize = false;
    bool preCompute = false;
    double thresh = 0.0;

    sgpp::datadriven::MultipleClassRefinementFunctor mcrf(grids2, alphas2,
            numRefinements, levelPenalize, preCompute, thresh);
    sgpp::datadriven::MultiGridRefinementFunctor* multifun = &mcrf;
    mcrf.printPointsPlott();
    std::cout << std::endl;

    std::cout << "Size: " << mcrf.getCombinedGrid()->getStorage().getSize() << std::endl;
    std::cout << "Dim: " << mcrf.getCombinedGrid()->getDimension() << std::endl;
    std::cout << std::endl;

    std::cout << "Evaluation:" << std::endl << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "   0   | " << eval.at(0) << " | " << eval.at(1) << std::endl;
    for ( size_t i = 1; i < numSteps + 1; i++ ) {
        mcrf.refine();

        // Re-train the grids. This has to happend AFTER all grids are refined
        for ( int i = 0 ; i < classes ; i++ ) {
            learner.at(i).train(*grids2.at(i), *alphas2.at(i), dataCl.at(i), lambda);
        }

        mcrf.prepareGrid(grids2, alphas2);
        std::cout << std::endl;
        eval = doClassification(grids2, alphas2, dataTrain, targetTrain);
        mcrf.printPointsPlott();
        mcrf.printPointsInfo();
        std::cout << "   " << i << "   | " << eval.at(0) << " | " << eval.at(1)
             << std::endl;
    }

    // double totalScore = multifun->getTotalRefinementValue(mcrf.getCombinedGrid()->getStorage());
    // std::cout << "TotalRefinementValue: " << totalScore << std::endl;

    return 0;
}

sgpp::datadriven::LearnerSGDE createSGDELearner(size_t dim, size_t level,
                                                double lambda) {
    sgpp::base::RegularGridConfiguration gridConfig;
    gridConfig.dim_ = dim;
    gridConfig.level_ = static_cast<int>(level);
    gridConfig.type_ = sgpp::base::GridType::Linear;

    // configure adaptive refinement
    sgpp::base::AdpativityConfiguration adaptConfig;
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
    regularizationConfig.regType_ =
        sgpp::datadriven::RegularizationType::Laplace;

    // configure learner
    sgpp::datadriven::CrossvalidationForRegularizationConfiguration
    crossvalidationConfig;
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

sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig,
                solverConfig, regularizationConfig, crossvalidationConfig);
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
      std::unique_ptr<sgpp::base::OperationEval>
        e(sgpp::op_factory::createOperationEval(*grids.at(i)));
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
      double ce = (100.0 * (1.0 - (static_cast<double>(classErrorCounts.at(i)) /
                                   classCounts.at(i))));
      ss << std::fixed << std::setprecision(2) << ce;
      if (i < grids.size() - 1) {
        ss << "  ";
      }
    }
    std::stringstream ss2;
    ss2 << std::fixed << std::setprecision(3);
    ss2 << (100.0 * (1.0 - (static_cast<double>(totalCount) /
                            static_cast<double>(testData.getNrows()))));
    std::vector<std::string> result;
    result.push_back(ss.str());
    result.push_back(ss2.str());
    return result;
}
