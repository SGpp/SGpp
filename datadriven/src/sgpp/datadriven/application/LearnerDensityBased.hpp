// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERDENSITYBASED_HPP_
#define LEARNERDENSITYBASED_HPP_

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/application/LearnerBase.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/globaldef.hpp>

#include <map>
#include <vector>

namespace sgpp {
namespace datadriven {

class LearnerDensityBased : public datadriven::LearnerBase {
 protected:
  // Mapping from class index to class number:
  std::map<int, double> indexToClass;
  // Stores the coefficients for every class
  std::vector<base::DataVector> alphaVec;
  /// regularization mode
  datadriven::RegularizationType CMode;
  // with prior
  bool withPrior;
  // number of classes
  size_t nrClasses;
  // prior of data
  std::vector<double> prior;
  // vectors of grids
  std::vector<std::unique_ptr<base::Grid>> gridVec;
  // vector of regterms
  std::vector<std::unique_ptr<base::OperationMatrix>> CVec;

 public:
  LearnerDensityBased(datadriven::RegularizationType&, const bool isRegression,
                      const bool isVerbose = true);
  virtual ~LearnerDensityBased();

  /**
   * Create a grid for each class
   *
   * @param GridConfig grid config
   */
  virtual void InitializeGrid(const base::RegularGridConfiguration& GridConfig);

  /**
   * Learning a dataset with spatially adaptive sparse grids
   *
   * @param testDataset the training dataset
   * @param classes classes corresponding to the training dataset
   * @param GridConfig configuration of the regular start grid
   * @param SolverConfigRefine configuration of the SLE solver during the adaptive refinements of
   *   the grid
   * @param SolverConfigFinal configuration of the final SLE solving step on the refined grid
   * @param AdaptConfig configuration of the adaptivity strategy
   * @param testAccDuringAdapt set to true if the training accuracy should be determined in evert
   *   refinement step
   * @param lambda regularization parameter lambda
   */
  virtual LearnerTiming train(base::DataMatrix& testDataset, base::DataVector& classes,
                              const base::RegularGridConfiguration& GridConfig,
                              const solver::SLESolverConfiguration& SolverConfigRefine,
                              const solver::SLESolverConfiguration& SolverConfigFinal,
                              const base::AdpativityConfiguration& AdaptConfig,
                              bool testAccDuringAdapt, const double lambda);

  virtual base::DataVector predict(base::DataMatrix& testDataset);
  /// construct system matrix
  std::unique_ptr<datadriven::DMSystemMatrixBase> createDMSystem(base::DataMatrix& trainDataset,
                                                                 double lambda) override;

  /**
   * Returns number of grid points for the density
   * with the maximum number of grid points
   */
  size_t getNrGridPoints();

  /**
   * Get Prior
   */
  bool getWithPrior() { return withPrior; }

  /**
   * Set prior
   *
   * @param p prior
   */
  bool setWithPrior(bool p) {
    withPrior = p;
    return withPrior;
  }

  /**
   * Get number of classes
   */
  size_t getNrClasses() { return nrClasses; }

  /**
   * Set number of classes
   *
   * @param c set number of classes
   */
  size_t setNrClasses(size_t c) {
    nrClasses = c;
    return nrClasses;
  }
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERDENSITYBASED_HPP_ */
