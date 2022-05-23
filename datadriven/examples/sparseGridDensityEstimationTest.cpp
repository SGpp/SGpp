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
#include <sgpp/datadriven/application/SparseGridDensityEstimator.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <random>
#include <string>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;

/**
 * First define an auxiliary function randu that returns a random point
 * (uniform, normal).
 */
void randu(DataVector& rvar, std::mt19937& generator) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    rvar[j] = distribution(generator);
  }
}

void randn(DataVector& rvar, std::mt19937& generator) {
  std::normal_distribution<double> distribution(0.5, 0.1);
  for (size_t j = 0; j < rvar.getSize(); ++j) {
    double value = -1.0;
    while (value <= 0.0 || value >= 1.0) {
      value = distribution(generator);
      rvar[j] = value;
    }
  }
}

/**
 * Second define another auxiliary function that calls the one defined above
 * multiple times and returns a matrix of random points (uniform, normal)
 */
void randu(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(static_cast<std::mt19937::result_type>(seedValue));
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randu(sample, generator);
    rvar.setRow(i, sample);
  }
}

void randn(DataMatrix& rvar, std::uint64_t seedValue = std::mt19937_64::default_seed) {
  size_t nsamples = rvar.getNrows(), ndim = rvar.getNcols();

  std::mt19937 generator(static_cast<std::mt19937::result_type>(seedValue));
  DataVector sample(ndim);
  for (size_t i = 0; i < nsamples; ++i) {
    randn(sample, generator);
    rvar.setRow(i, sample);
  }
}

/**
 * Now the main function begins by loading the test data from a file specified
 * in the string filename.
 */
int main(int argc, char** argv) {
  /**
   * Define number of dimensions of the toy problem.
   */
  size_t numDims = 3;

  /**
   * Load normally distributed samples.
   */
  sgpp::base::DataMatrix samples(1000, numDims);
  randn(samples);

  /**
   * Configure the sparse grid of level 3 with linear basis functions and the
   * same dimension as the given test data. \n
   * Alternatively load a sparse grid that has been saved to a file, see the
   * commented line.
   */
  std::cout << "# create grid config" << std::endl;
  sgpp::base::RegularGridConfiguration gridConfig;
  gridConfig.dim_ = numDims;
  gridConfig.level_ = 3;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  /**
   * Configure the adaptive refinement. Therefore the number of refinements and
   * the number of points are specified.
   */
  std::cout << "# create adaptive refinement config" << std::endl;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.numRefinementPoints_ = 10;

  /**
   * Configure the solver. The solver type is set to the conjugent gradient
   * method and the maximum number of iterations, the tolerance epsilon and the
   * threshold are specified.
   */
  std::cout << "# create solver config" << std::endl;
  sgpp::solver::SLESolverConfiguration solverConfig;
  solverConfig.type_ = sgpp::solver::SLESolverType::CG;
  solverConfig.maxIterations_ = 1000;
  solverConfig.eps_ = 1e-14;
  solverConfig.threshold_ = 1e-14;
  solverConfig.verbose_ = true;

  /**
   * Configure the regularization for the laplacian operator.
   */
  std::cout << "# create regularization config" << std::endl;
  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.type_ = sgpp::datadriven::RegularizationType::Laplace;

  /**
   * Configure the learner by specifying: \n
   * - enable kfold, \n
   * - an initial value for the Lagrangian multiplier \f$\lambda\f$ and the
   * interval \f$ [\lambda_{Start} , \lambda_{End}] \f$ in which \f$\lambda\f$
   * will be searched, \n
   * - whether a logarithmic scale is used, \n
   * - the parameters shuffle and an initial seed for the random value
   * generation, \n
   * - whether parts of the output shall be kept off.
   */
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

  // configure learner
  std::cout << "# create learner config" << std::endl;
  sgpp::datadriven::SGDEConfiguration sgdeConfig;
  sgdeConfig.makePositive_ = true;
  sgdeConfig.makePositive_candidateSearchAlgorithm_ =
      sgpp::datadriven::MakePositiveCandidateSearchAlgorithm::HybridFullIntersections;
  sgdeConfig.makePositive_interpolationAlgorithm_ =
      sgpp::datadriven::MakePositiveInterpolationAlgorithm::InterpolateBoundaries1d;
  sgdeConfig.unitIntegrand_ = true;

  /**
   * Create the learner using the configurations set above. Then initialize it
   * with the data read from the file in the first step and train the learner.
   */
  std::cout << "# creating the learner" << std::endl;
  sgpp::datadriven::SparseGridDensityEstimator learner(gridConfig, adaptivityConfig, solverConfig,
                                                       regularizationConfig, crossvalidationConfig,
                                                       sgdeConfig);
  learner.initialize(samples);

  /**
   * Estimate the probability density function (pdf) via a Gaussian kernel
   * density estimation (KDE) and print the corresponding values.
   */
  sgpp::datadriven::KernelDensityEstimator kde(samples);
  sgpp::base::DataVector x(learner.getDim(), 0.5);

  std::cout << "--------------------------------------------------------\n";
  std::cout << learner.getSurpluses().getSize() << " -> " << learner.getSurpluses().sum() << "\n";
  std::cout << "pdf_SGDE(x) = " << learner.pdf(x) << " ~ " << kde.pdf(x) << " =pdf_KDE(x)\n";
  std::cout << "mean_SGDE(x) = " << learner.mean() << " ~ " << kde.mean() << " = mean_KDE(x)\n";
  std::cout << "var_SGDE(x) = " << learner.variance() << " ~ " << kde.variance() << "=var_KDE(x)\n";

  sgpp::base::DataMatrix* bounds = new DataMatrix(gridConfig.dim_, 2);
  for (size_t idim = 0; idim < gridConfig.dim_; idim++) {
    bounds->set(idim, 0, 0.0);
    bounds->set(idim, 1, 1.0);
  }

  /**
   * Print the covariances.
   */
  sgpp::base::DataMatrix C(gridConfig.dim_, gridConfig.dim_);
  std::cout << "---------------------- Cov_SGDE ------------------------------" << std::endl;
  learner.cov(C, bounds);
  std::cout << C.toString() << std::endl;

  std::cout << "---------------------- Cov KDE--------------------------------" << std::endl;
  kde.cov(C);
  std::cout << C.toString() << std::endl;

  /**
   * Apply the inverse Rosenblatt transformation to a matrix of random points.
   * To do this first generate the random points via randu, then initialize an
   * inverse Rosenblatt transformation operation and apply it to the points.
   * Finally print the calculated values.
   */
  std::cout << "------------------------------------------------------" << std::endl;
  // inverse Rosenblatt transformation
  sgpp::datadriven::OperationInverseRosenblattTransformation* opInvRos(
      sgpp::op_factory::createOperationInverseRosenblattTransformation(learner.getGrid()));
  sgpp::base::DataMatrix points(12, gridConfig.dim_);
  randu(points);

  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "uniform space" << std::endl;
  std::cout << points.toString() << std::endl;

  sgpp::base::DataMatrix pointsCdf(points.getNrows(), points.getNcols());
  opInvRos->doTransformation(&learner.getSurpluses(), &points, &pointsCdf);

  /**
   * To check whether the results are correct, perform a Rosenblatt
   * transformation on the data that has been created by the inverse Rosenblatt
   * transformation above and print the calculated values.
   */
  points.setAll(0.0);
  sgpp::datadriven::OperationRosenblattTransformation* opRos(
      sgpp::op_factory::createOperationRosenblattTransformation(learner.getGrid()));
  opRos->doTransformation(&learner.getSurpluses(), &pointsCdf, &points);
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "original space" << std::endl;
  std::cout << pointsCdf.toString() << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "uniform space" << std::endl;
  std::cout << points.toString() << std::endl;
}
