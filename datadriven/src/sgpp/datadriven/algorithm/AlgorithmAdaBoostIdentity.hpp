// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ALGORITHMADABOOSTIDENTITY_HPP
#define ALGORITHMADABOOSTIDENTITY_HPP

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/*
 * Algorithm for Adaboosting
 *
 * For training the grid, a sparse grid standard approach is taken.
 */
class AlgorithmAdaBoostIdentity : public AlgorithmAdaBoostBase {
 protected:
  /**
   * Performs a solver to get alpha use DMWeightMatrix as the System Matrix
   *
   * @param lambda the regularisation parameter
   * @param weight the weights of examples
   * @param alpha output the coefficients of the sparse grid's basis functions
   * @param final judgement the final step of this base learner
   */
  virtual void alphaSolver(double& lambda, base::DataVector& weight, base::DataVector& alpha,
                           bool final);

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param gridType reference to the type of grid
   * @param gridLevel reference to the level of grid
   * @param trainData reference to the training dataset
   * @param trainDataClass reference to the class of training dataset
   * @param NUM the number of baselearner for Adaboosting
   * @param lambda the regularisation parameter
   * @param IMAX the parameter for ConjugateGradients
   * @param eps the parameter for ConjugateGradients
   * @param IMAX_final the parameter for ConjugateGradients used for last refinement step
   * @param eps_final the parameter for ConjugateGradients used for last refinement step
   * @param firstLabel one label from training dataset
   * @param secondLabel another label from training dataset
   * @param threshold the parameter for predicting a class
   * @param maxLambda the max lambda used in searching optimal lambda
   * @param minLambda the min lambda used in searching optimal lambda
   * @param searchNum the searching times used in searching for optimal lambda
   * @param refine the judgement of refine
   * @param refineMode Select the refine mode
   * @param refineNum the Number of refinement with a certain percentage of Grid points
   * @param numberOfAda the number of Grid points to refine
   * @param percentOfAda the percentage of Grid points to refine
   * @param mode the adaboost type to choose
   */
  AlgorithmAdaBoostIdentity(base::Grid& SparseGrid, base::GridType gridType,
                            base::level_t gridLevel, base::DataMatrix& trainData,
                            base::DataVector& trainDataClass, size_t NUM, double lambda,
                            size_t IMAX, double eps, size_t IMAX_final, double eps_final,
                            double firstLabel, double secondLabel, double threshold,
                            double maxLambda, double minLambda, size_t searchNum, bool refine,
                            size_t refineMode, size_t refineNum, size_t numberOfAda,
                            double percentOfAda, size_t mode);

  /**
   * Std-Deconstructor
   */
  virtual ~AlgorithmAdaBoostIdentity();
};

}  // namespace datadriven
}  // namespace sgpp
#endif
