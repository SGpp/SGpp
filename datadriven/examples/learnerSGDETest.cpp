// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <string>
#include <random>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/operation/hash/OperationMakePositive.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/datadriven/application/LearnerSGDE.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/datadriven/application/RegularizationConfiguration.hpp"
#include "sgpp/datadriven/application/KernelDensityEstimator.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;

void randu(DataVector& rvar, std::mt19937& generator) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
  }
}

void randu(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(seedValue);
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randu(sample, generator);
    rvar.setRow(i, sample);
  }
}

int main(int argc, char** argv) {
  std::string filename = "../../parallel/tests/data/friedman_4d_2000.arff";

  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset dataset = sgpp::datadriven::ARFFTools::readARFF(filename);
  sgpp::base::DataMatrix& samples = dataset.getData();

  // configure grid
  std::cout << "# create grid config" << std::endl;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = dataset.getDimension();
  gridConfig.level_ = 2;
  gridConfig.type_ = sgpp::base::GridType::Linear;
  gridConfig.maxDegree_ = 3;
  //  gridConfig.filename_ = "/tmp/sgde-grid-4391dc6e-54cd-4ca2-9510-a9c02a2889ec.grid";

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
  solverConfig.eps_ = 1e-14;
  solverConfig.threshold_ = 1e-14;
  solverConfig.verbose_ = true;

  // configure regularization
  std::cout << "# create regularization config" << std::endl;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.regType_ = sgpp::datadriven::RegularizationType::Laplace;

  // configure cross validation
  std::cout << "# create cross validation config" << std::endl;
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

  // configure learner
  std::cout << "# create learner config" << std::endl;
  sgpp::datadriven::SGDEConfiguration sgdeConfig;
  sgdeConfig.makePositive_ = false;
  sgdeConfig.makePositive_candidateSearchAlgorithm_ =
      sgpp::base::MakePositiveCandidateSearchAlgorithm::Intersections;
  sgdeConfig.makePositive_interpolationAlgorithm_ =
      sgpp::base::MakePositiveInterpolationAlgorithm::SetToZero;
  sgdeConfig.makePositive_verbose_ = false;
  sgdeConfig.unitIntegrand_ = false;

  std::cout << "# creating the learner" << std::endl;
  sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
                                        crossvalidationConfig, sgdeConfig);
  learner.initialize(samples);

  std::cout << "# estimating a kde" << std::endl;
  sgpp::datadriven::KernelDensityEstimator kde(samples);

  sgpp::base::DataVector x(learner.getDim());

  for (size_t i = 0; i < x.getSize(); i++) {
    x[i] = 0.5;
  }

  std::cout << "--------------------------------------------------------" << std::endl;
  std::cout << learner.getSurpluses().getSize() << " -> " << learner.getSurpluses().sum()
            << std::endl;
  std::cout << "pdf_SGDE(x) = " << learner.pdf(x) << " ~ " << kde.pdf(x) << " = pdf_KDE(x)"
            << std::endl;
  /*std::cout << "mean_SGDE(x) = " << learner.mean() << " ~ " << kde.mean() << " = mean_KDE(x)"
            << std::endl;
  std::cout << "var_SGDE(x) = " << learner.variance() << " ~ " << kde.variance() << " = var_KDE(x)"
            << std::endl;*/

  /*sgpp::base::DataMatrix C(gridConfig.dim_, gridConfig.dim_);
  std::cout << "---------------------- Cov_SGDE ------------------------------" << std::endl;
  learner.cov(C);
  std::cout << C.toString() << std::endl;

  std::cout << "---------------------- Cov KDE--------------------------------" << std::endl;
  kde.cov(C);
  std::cout << C.toString() << std::endl;

  std::cout << "------------------------------------------------------" << std::endl;*/

  // inverse Rosenblatt transformation
  auto opInvRos =
      sgpp::op_factory::createOperationInverseRosenblattTransformation(learner.getGrid());
  sgpp::base::DataMatrix points(12, gridConfig.dim_);
  randu(points);

  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << points.toString() << std::endl;

  sgpp::base::DataMatrix pointsCdf(points.getNrows(), points.getNcols());
  opInvRos->doTransformation(&learner.getSurpluses(), &points, &pointsCdf);

  points.setAll(0.0);
  auto opRos = sgpp::op_factory::createOperationRosenblattTransformation(learner.getGrid());
  opRos->doTransformation(&learner.getSurpluses(), &pointsCdf, &points);
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << pointsCdf.toString() << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << points.toString() << std::endl;
}
