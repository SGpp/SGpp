// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/MultipleClassPoint.hpp>

#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/MultipleClassRefinementFunctor.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <sys/resource.h>
#include <cmath>
#include <ctime>

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <locale>
#include <chrono>

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
                                          sgpp::base::DataVector& testLabel,
                                          size_t classes);
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


int main(int argc, char* argv[]) {
    bool useZCRF = false;
    size_t x = 4;
    size_t y = 2;
    size_t z = 2;
    if ( argc == 5 ) {
        useZCRF = std::atoi(argv[1]) == 1;
        x = std::atoi(argv[2]);
        y = std::atoi(argv[3]);
        z = std::atoi(argv[4]);
    }
/*
    std::string filename2 = "exDim4Class" + std::to_string(z);
    size_t classes = z;
    size_t dim = 4;
    size_t level = 4;
    double lambda = 1e-2;
    size_t numSteps = 8;
    size_t numRefinements = x;
    size_t partCombined = y;
//*/
///*
    std::string filename2 = "mulitpleClassesTest";
    size_t classes = 4;
    size_t dim = 2;
    size_t level = 4;
    double lambda = 1e-2;
    size_t numSteps = 5;
    size_t numRefinements = 3;
    size_t partCombined = 0;
//*/
/*
    std::string filename2 = "mulitpleClassesTest2";
    size_t classes = 4;
    size_t dim = 2;
    size_t level = 3;
    double lambda = 1e-2;
    size_t numSteps = 6;
    size_t numRefinements = 4;
    size_t partCombined = 1;
//*/
/*
    std::string filename2 = "starsgalaxies";
    size_t classes = 3;
    size_t dim = 4;
    size_t level = 5;
    double lambda = 1e-7;
    size_t numSteps = 8;
    size_t numRefinements = 4;
    size_t partCombined = 2;
//*/
/*
    std::string filename2 = "star_test";
    size_t classes = 2;
    size_t dim = 4;
    size_t level = 3;
    double lambda = 1e-5;
    size_t numSteps = 8;
    size_t numRefinements = 3;
    size_t partCombined = 0;
//*/
    if ( argc == 5 ) {
        numRefinements = x;
        partCombined = y;
    }
    //------------- only calculation after -------------------------
    std::time_t starttime = std::time(nullptr);

   /* sgpp::datadriven::Dataset dataset =
            sgpp::datadriven::ARFFTools::readARFF(
            "../../datadriven/tests/data/" + filename2 + ".arff");*/
    sgpp::datadriven::Dataset dataset =
            sgpp::datadriven::ARFFTools::readARFF(
            "/home/katrin/Desktop/multiClass/" + filename2 + ".arff");
    sgpp::base::DataMatrix dataTrain = dataset.getData();
    sgpp::base::DataVector targetTrain = dataset.getTargets();

    // write log to file
    std::string filename = "/home/katrin/Desktop/multi_log_"
            + filename2 + "_" + std::to_string(numRefinements)
            + "_" + std::to_string(partCombined)
            + "_" + std::to_string(useZCRF) +  ".log";
    std::ofstream out(filename);
    std::cout.rdbuf(out.rdbuf());

    std::cout << "Read training data: " << dataTrain.getNrows() << std::endl;

    std::vector<sgpp::base::DataMatrix> dataCl;
    std::vector<sgpp::datadriven::LearnerSGDE> learner;
    for ( size_t i = 0 ; i < classes ; i++ ) {
        dataCl.push_back(sgpp::base::DataMatrix(0.0, dataTrain.getNcols()));
    }

    sgpp::base::DataVector row(dataTrain.getNcols());
    for ( size_t i = 0 ; i < dataTrain.getNrows() ; i++ ) {
        dataTrain.getRow(i, row);
        dataCl.at((size_t)targetTrain.get(i)).appendRow(row);
    }

    /**
    * Approximate a probability density function for the class data using
    * LearnerSGDE, one for each class. Initialize the learners with the data
    */
    for ( size_t i = 0 ; i < classes ; i++ ) {
        std::cout << "Data points of class " << std::setw(3) << std::right << i << ": ";
        std::cout << std::setw(14) << std::right << dataCl.at(i).getNrows() << "   | ";
        learner.push_back(createSGDELearner(dim, level, lambda));
        learner.back().initialize(dataCl.at(i));
    }

    /**
     * Bundle grids and surplus vector pointer needed for refinement and
     * evaluation
     */
    std::vector<sgpp::base::Grid*> grids2;
    std::vector<sgpp::base::DataVector*> alphas2;
    for ( size_t i = 0 ; i < classes ; i++ ) {
        grids2.push_back(learner.at(i).getGrid().get());
        alphas2.push_back(learner.at(i).getSurpluses().get());
    }

    std::time_t midtime = std::time(nullptr);
    std::cout << "needed time to read data: " << midtime - starttime  << "s" << std::endl;
    std::cout << "Start refining" << std::endl;

    std::cout << "DATA SETTINGS:" << std::endl;
    std::cout << "  filename:    " << filename2 << std::endl;
    std::cout << "  classes:     " << classes << std::endl;
    std::cout << "  dimensions:  " << dim << std::endl;
    std::cout << "  start level: " << level << std::endl;
    std::cout << "  lambda:      " << lambda << std::endl;
    std::cout << "  numRef:      " << numRefinements << std::endl;
    std::cout << "    partComb:  " << partCombined << std::endl;
    if ( useZCRF ) {
        std::cout << "Using the ZeroCrossingRefinementFunctor" << std::endl;
    } else {
        std::cout << "Using the MultipleClassRefinementFunctor" << std::endl;
    }
    std::cout << "---------------------------------------------" << std::endl;

    auto c_start = std::chrono::high_resolution_clock::now();
    auto c_ref = std::chrono::high_resolution_clock::now();
    std::cout << "Sizes:";
    for (size_t i = 0; i < classes ; i++) {
        std::cout << " (" << i << "->" << grids2.at(i)->getStorage().getSize() << ") ";
    }
    std::cout << std::endl;
    auto c_learn = std::chrono::high_resolution_clock::now();
    for ( size_t i = 0 ; i < classes ; i++ ) {
        learner.at(i).train(*grids2.at(i), *alphas2.at(i), dataCl.at(i), lambda);
    }
    auto c_eval = std::chrono::high_resolution_clock::now();
    std::vector<std::string> eval = doClassification(grids2, alphas2,
                    dataTrain, targetTrain, classes);
    std::cout << "   0   | " << eval.at(0) << " | " << eval.at(1) << std::endl;
    auto c_end = std::chrono::high_resolution_clock::now();
    std::cout << "time: 0: sum: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(c_end-c_start).count()
          << " ms  | step: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(c_end-c_ref).count()
          << " ms  | ref: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(c_learn-c_ref).count()
          << " ms  | learn: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(c_eval-c_learn).count()
          << " ms  | eval: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(c_end-c_eval).count()
          << " ms" << std::endl;

    // amount of points refined in the single classes vs combined grid
    // class grid: numRefinements - partCombined
    // combined grid: partCombined
    bool levelPenalize = true;
    bool preCompute = false;
    double thresh = 0;

    // Zero-crossing-based refinement
    sgpp::datadriven::MultiGridRefinementFunctor* fun = nullptr;
    sgpp::datadriven::ZeroCrossingRefinementFunctor funZrcr(grids2, alphas2,
            numRefinements, levelPenalize, preCompute);
    fun = &funZrcr;

    sgpp::datadriven::MultipleClassRefinementFunctor* multifun = nullptr;
    sgpp::datadriven::MultipleClassRefinementFunctor mcrf(grids2, alphas2,
            numRefinements, thresh);
    multifun = &mcrf;
    multifun->printPointsPlott();

    for ( size_t i = 1; i < numSteps + 1; i++ ) {
        std::cout << "---------------------------------------------" << std::endl;
        c_ref = std::chrono::high_resolution_clock::now();
        if ( useZCRF ) {
            for ( size_t i = 0 ; i < classes ; i++ ) {
                fun->setGridIndex(i);
                grids2.at(i)->getGenerator().refine(*fun);
            }
        } else {
            multifun->refine(partCombined);
        }
        std::cout << "Sizes:";
        for (size_t i = 0; i < classes ; i++) {
            std::cout << " (" << i << "->" << grids2.at(i)->getStorage().getSize() << ") ";
        }
        std::cout << std::endl;

        c_learn = std::chrono::high_resolution_clock::now();
        // Re-train the grids. This has to happend AFTER all grids are refined
        for ( size_t i = 0 ; i < classes ; i++ ) {
            learner.at(i).train(*grids2.at(i), *alphas2.at(i), dataCl.at(i), lambda);
        }

        c_eval = std::chrono::high_resolution_clock::now();
        eval = doClassification(grids2, alphas2, dataTrain, targetTrain, classes);
        std::cout << "   " << i << "   | " << eval.at(0) << " | " << eval.at(1)
             << std::endl;

        c_end = std::chrono::high_resolution_clock::now();
        std::cout << "time: " << i << ": sum: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(c_end-c_start).count()
              << " ms  | step: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(c_end-c_ref).count()
              << " ms  | ref: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(c_learn-c_ref).count()
              << " ms  | learn: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(c_eval-c_learn).count()
              << " ms  | eval: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(c_end-c_eval).count()
              << " ms" << std::endl;
    }

    // if (dim == 2) {
    //     multifun->printPointsPlott();
    // }
    // multifun->printPointsPlott();
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "rusage: memory: " << usage.ru_maxrss
        << " kB  | time user: " << usage.ru_utime.tv_sec
        << " s  | time system: " << usage.ru_stime.tv_sec
        << " s" << std::endl;
    multifun->printPointsInfo();
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
                  sgpp::base::DataVector& testLabel,
                  size_t classes) {
    double best_eval = -1000.0;
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
    std::vector<std::vector<int> > gridEval(classes+1, std::vector<int>(classes+1, 0));
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
      // set array
      gridEval[static_cast<size_t>(floor(testLabel.get(i)))]
                    [static_cast<size_t>(indices.get(i))] += 1;
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

    bool printMapping = false;
    if ( printMapping ) {
        // print array
        std::cout << std::setw(10) << std:: left << "classes" << " | ";
        std::cout << std::setw(10) << std:: left << "all Points" << " | ";
        for (size_t a = 0; a < grids.size(); a++) {
          std::cout << "map in " << std::setw(3) << std:: left << a << " | ";
        }
        std::cout << std::setw(10) << std:: left << "wrong map" << std::endl;
        for (size_t a = 0; a < grids.size(); a++) {
          std::cout << " Class: " << std::setw(2) << std:: left << a << " | ";
          std::cout << std::setw(10) << std::right << classCounts.at(a) << " | ";
          for (size_t b = 0; b < grids.size(); b++) {
            std::cout << std::setw(10) << std::right << gridEval[a][b] << " | ";
          }
          std::cout << std::setw(10) << std::right << classErrorCounts.at(a) << std::endl;
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
