// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_learnerSGDETest_cpp learner SGDE
 * This examples demonstrates density estimation.
 */

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/application/KernelDensityEstimator.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>

#include <random>
#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;

/**
 * First define an auxiliary function randu that returns a random point.
 */
void randu(DataVector& rvar, std::mt19937& generator) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
  }
}

/**
 * Second define another auxiliary function that calls the one defined above multiple times and
 * returns a matrix of random points.
 */
void randu(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(seedValue);
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randu(sample, generator);
    rvar.setRow(i, sample);
  }
}

/**
 * Now the main function begins by loading the test data from a file specified in the string
 * filename.
 */
int main(int argc, char** argv) {
  std::string filename = "../datasets/friedman/friedman2_4d_10000.arff";

  std::cout << "# loading file: " << filename << std::endl;
  sgpp::datadriven::Dataset dataset =
    sgpp::datadriven::ARFFTools::readARFFFromFile(filename);
  sgpp::base::DataMatrix& samples = dataset.getData();

  /**
   * Configure the sparse grid of level 3 with linear basis functions and the same dimension as the
   * given test data.
   * Alternatively load a sparse grid that has been saved to a file, see the commented line.
   */
  std::cout << "# create grid config" << std::endl;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = dataset.getDimension();
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;
  //  gridConfig.filename_ = "/tmp/sgde-grid-4391dc6e-54cd-4ca2-9510-a9c02a2889ec.grid";

  /**
   * Configure the adaptive refinement. Therefore the number of refinements and the number of points
   * are specified.
   */
  std::cout << "# create adaptive refinement config" << std::endl;
  sgpp::base::AdaptivityConfiguration adaptConfig;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.noPoints_ = 10;

  /**
   * Configure the solver. The solver type is set to the conjugent gradient method and the maximum
   * number of iterations, the tolerance epsilon and the threshold are specified.
   */
  std::cout << "# create solver config" << std::endl;
  sgpp::solver::SLESolverConfiguration solverConfig;
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-14;
  solverConfig.threshold_ = 1e-14;

  /**
   * Configure the regularization for the laplacian operator.
   */
  std::cout << "# create regularization config" << std::endl;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Laplace;

  /**
   * Configure the learner by specifying:
   * - an initial value for the lagrangian multiplier \f$\lambda\f$ and the interval
   *   \f$ [\lambda_{Start} , \lambda_{End}] \f$ in which \f$\lambda\f$ will be searched,
   * - whether a logarithmic scale is used,
   * - the parameters shuffle and an initial seed for the random value generation,
   * - whether parts of the output shall be kept off.
   */
  std::cout << "# create learner config" << std::endl;
  sgpp::datadriven::CrossvalidationConfiguration crossvalidationConfig;
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

  /**
   * Create the learner using the configuratons set above. Then initialize it with the data read
   * from the file in the first step and train the learner.
   */
  std::cout << "# creating the learner" << std::endl;
  sgpp::datadriven::LearnerSGDE learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
                                        crossvalidationConfig);
  learner.initialize(samples);
  learner.train();
  /**
   * Estimate the probability density function (pdf) via a gaussian kernel density estimation (KDE)
   * and print the corresponding values.
   */
  sgpp::datadriven::KernelDensityEstimator kde(samples);
  sgpp::base::DataVector x(learner.getDim());
  x.setAll(0.5);

  std::cout << "--------------------------------------------------------\n";
  std::cout << learner.getSurpluses()->getSize() << " -> " << learner.getSurpluses()->sum() << "\n";
  std::cout << "pdf_SGDE(x) = " << learner.pdf(x) << " ~ " << kde.pdf(x) << " =pdf_KDE(x)\n";
  std::cout << "mean_SGDE(x) = " << learner.mean() << " ~ " << kde.mean() << " = mean_KDE(x)\n";
  std::cout << "var_SGDE(x) = " << learner.variance() << " ~ " << kde.variance() << "=var_KDE(x)\n";

  /**
   * Print the covariances.
   */
  sgpp::base::DataMatrix C(gridConfig.dim_, gridConfig.dim_);
  std::cout << "---------------------- Cov_SGDE ------------------------------" << std::endl;
  learner.cov(C);
  std::cout << C.toString() << std::endl;

  std::cout << "---------------------- Cov KDE--------------------------------" << std::endl;
  kde.cov(C);
  std::cout << C.toString() << std::endl;

  /**
   * Apply the inverse Rosenblatt transformation to a matrix of random points. To do this first
   * generate the random points via randu, then initialize an inverse Rosenblatt transformation
   * operation and apply it to the points. Finally print the calculated values.
   */
  std::cout << "------------------------------------------------------" << std::endl;
  // inverse Rosenblatt transformation
  sgpp::datadriven::OperationInverseRosenblattTransformation* opInvRos(
      sgpp::op_factory::createOperationInverseRosenblattTransformation(*learner.getGrid()));
  sgpp::base::DataMatrix points(12, gridConfig.dim_);
  randu(points);

  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << points.toString() << std::endl;

  sgpp::base::DataMatrix pointsCdf(points.getNrows(), points.getNcols());
  opInvRos->doTransformation(learner.getSurpluses(), &points, &pointsCdf);

  /**
   * To check whether the results are correct perform a Rosenform transformation on the data that
   * has been created by the inverse Rosenblatt transformation above and print the calculated
   * values.
   */
  points.setAll(0.0);
  sgpp::datadriven::OperationRosenblattTransformation* opRos(
      sgpp::op_factory::createOperationRosenblattTransformation(*learner.getGrid()));
  opRos->doTransformation(learner.getSurpluses(), &pointsCdf, &points);
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << pointsCdf.toString() << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << points.toString() << std::endl;
}
