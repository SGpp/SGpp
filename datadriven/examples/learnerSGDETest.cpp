// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <string>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/application/LearnerSGDE.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/datadriven/application/RegularizationConfiguration.hpp"
#include "sgpp/datadriven/application/GaussianKDE.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"

int main(int argc, char** argv) {
  std::string filename = "../tests/data/friedman_4d_2000.arff";

  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset dataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& samples = dataset.getData();

  // configure grid
  std::cout << "# create grid config" << std::endl;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = dataset.getDimension();
  gridConfig.level_ = 2;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // configure adaptive refinement
  std::cout << "# create adaptive refinement config" << std::endl;
  sgpp::base::AdpativityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.noPoints_ = 10;

  // configure solver
  std::cout << "# create solver config" << std::endl;
  sgpp::solver::SLESolverConfiguration solverConfig;
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-10;
  solverConfig.threshold_ = 1e-10;

  // configure regularization
  std::cout << "# create regularization config" << std::endl;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;

  // configure learner
  std::cout << "# create learner config" << std::endl;
  sgpp::datadriven::CrossvalidationForRegularizationConfiguration crossvalidationConfig;
  crossvalidationConfig.enable_ = false;
  crossvalidationConfig.kfold_ = 3;
  crossvalidationConfig.lambda_ = 3.16228e-06;
  crossvalidationConfig.lambdaStart_ = 1e-1;
  crossvalidationConfig.lambdaEnd_ = 1e-10;
  crossvalidationConfig.lambdaSteps_ = 3;
  crossvalidationConfig.logScale_ = true;
  crossvalidationConfig.shuffle_ = true;
  crossvalidationConfig.seed_ = 1234567;
  crossvalidationConfig.silent_ = false;

  std::cout << "# creating the learner" << std::endl;
  sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
                                        crossvalidationConfig);
  learner.initialize(samples);

  sgpp::base::DataMatrix covSgde(learner.getDim(), learner.getDim());
  learner.cov(covSgde);
  std::cout << covSgde.toString() << std::endl;

  sgpp::datadriven::GaussianKDE kde(samples);
  sgpp::base::DataMatrix covKDE(kde.getDim(), kde.getDim());
  kde.cov(covKDE);
  std::cout << covKDE.toString() << std::endl;

  sgpp::base::DataVector x(learner.getDim());

  for (size_t i = 0; i < x.getSize(); i++) {
    x[i] = 0.5;
  }

  std::cout << "--------------------------------------------------------" << std::endl;
  std::cout << learner.getSurpluses()->getSize() << " -> " << learner.getSurpluses()->sum()
            << std::endl;
  std::cout << "pdf_SGDE(x) = " << learner.pdf(x) << " ~ " << kde.pdf(x) << " = pdf_KDE(x)"
            << std::endl;
  std::cout << "mean_SGDE(x) = " << learner.mean() << " ~ " << kde.mean() << " = mean_KDE(x)"
            << std::endl;
  std::cout << "var_SGDE(x) = " << learner.variance() << " ~ " << kde.variance() << " = var_KDE(x)"
            << std::endl;

  sgpp::base::DataMatrix C(gridConfig.dim_, gridConfig.dim_);
  std::cout << "---------------------- Cov_SGDE ------------------------------" << std::endl;
  learner.cov(C);
  std::cout << C.toString() << std::endl;

  std::cout << "---------------------- Cov KDE--------------------------------" << std::endl;
  kde.cov(C);
  std::cout << C.toString() << std::endl;
}
