// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERDENSITYBASEDREG_HPP_
#define LEARNERDENSITYBASEDREG_HPP_

#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Class that implements a regression learner
 * using density estimation on Sparse Grids.
 *
 */
class LearnerDensityBasedReg: public LearnerBase {
 public:
  using LearnerBase::predict;

  /**
   * Constructor
   *
   * @param regularization regularization type
   * @param border offset for the normalization of the class data
   */
  LearnerDensityBasedReg(
    sgpp::datadriven::RegularizationType& regularization,
    double border = 0.);

  /**
   * Destructor
   */
  virtual ~LearnerDensityBasedReg();

  /**
   * Learning a dataset with spatially adaptive sparse grids
   *
   * @param trainDataset the training dataset
   * @param classes classes corresponding to the training dataset
   * @param GridConfig configuration of the regular start grid
   * @param SolverConfigRefine configuration of the SLE solver during the adaptive refinements of the grid
   * @param SolverConfigFinal configuration of the final SLE solving step on the refined grid
   * @param AdaptConfig configuration of the adaptivity strategy
   * @param testAccDuringAdapt set to true if the training accuracy should be determined in evert refinement step
   * @param lambda regularization parameter lambda
   */
  virtual LearnerTiming train(sgpp::base::DataMatrix& trainDataset,
                              sgpp::base::DataVector& classes,
                              const sgpp::base::RegularGridConfiguration& GridConfig,
                              const sgpp::solver::SLESolverConfiguration& SolverConfigRefine,
                              const sgpp::solver::SLESolverConfiguration& SolverConfigFinal,
                              const sgpp::base::AdpativityConfiguration& AdaptConfig,
                              bool testAccDuringAdapt, const double lambda);

  /**
   * Executes a regression test for a given dataset and returns the result
   *
   * @param testDataset dataset that is evaluated with the current learner
   *
   * @return regression values of testDataset
   */
  virtual sgpp::base::DataVector predict(sgpp::base::DataMatrix& testDataset);

  /**
   * simple dump of the one dimensional density function for a specific data point
   *  into file, e.g. used to plot with gnuplot.
   *
   * @param point point for which the one dimensional density function is computed
   * @param fileName filename to store the dump to
   * @param resolution resolution of function plot
   */
  void dumpDensityAtPoint(sgpp::base::DataVector& point, std::string fileName,
                          unsigned int resolution);

 protected:
  /// regularization mode
  sgpp::datadriven::RegularizationType CMode_;
  /// regularization operator
  sgpp::base::OperationMatrix* C_;
  /// maximum value (used for de-normalization)
  double maxValue_;
  /// minimum value (used for de-normalization)
  double minValue_;
  /// border for normalization of the class vector
  double border_;

  /**
   * inherited from LearnerBase, but not used
   */
  virtual sgpp::datadriven::DMSystemMatrixBase* createDMSystem(
    sgpp::base::DataMatrix& trainDataset, double lambda);
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERDENSITYBASEDREG_HPP_ */

