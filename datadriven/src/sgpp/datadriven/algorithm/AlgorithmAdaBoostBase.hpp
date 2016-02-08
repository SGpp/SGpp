// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ALGORITHMADABOOSTBASE_HPP
#define ALGORITHMADABOOSTBASE_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisation.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>

#include <cmath>
#include <string>
#include <utility>
#include <iostream>
#include <cstdlib>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

/*
 * Algorithm for Adaboosting
 * This algorithm is to train Train base learner according to sample distribution and obtain hypothesis
 * get the hypothesis weight
 * Then combine hypothesis linearly
 *
 * The main idea behind the algorithm is to get a better accuracy in classify dateset according to the training dataset
 *
 */
class AlgorithmAdaBoostBase {
 protected:
  /// the lambda, the regularisation parameter
  float_t lamb;
  /// the size of the grid
  size_t numData;
  /// Pointer to the data matrix
  SGPP::base::DataMatrix* data;
  /// DataSet Dimension
  size_t dim;
  /// Pointer to the class(real value in regression) of the data vector
  SGPP::base::DataVector* classes;
  /// Number of base learner for Adaboosting
  size_t numBaseLearners;
  /// the grid
  SGPP::base::Grid* grid;
  /// type of grid possible value are 1, 2 or 3(1 = Linear Grid, 2 = LinearL0Boundary Grid, 3 = ModLinear Grid);
  size_t type;
  /// Number of grid points
  size_t gridPoint;
  /// Number of the maximum grid points used in the algorithm
  SGPP::base::DataVector* maxGridPoint;
  /// Number of the sum grid points used in the algorithm
  SGPP::base::DataVector* sumGridPoint;
  /// grid level
  SGPP::base::HashGenerator::level_t level;
  /// Parameter for CG solver(during the refinement)
  size_t imax;
  /// Parameter for CG solver(during the refinement)
  float_t epsilon;
  /// Parameter for CG solver(for the last refinement)
  size_t imax_final;
  /// Parameter for CG solver(for the last refinement)
  float_t epsilon_final;
  /// One label of the DataSet
  float_t labelOne;
  /// Another label of the DataSet
  float_t labelTwo;
  /// Threshold to predict class
  float_t threshold;
  /// Log of the Max lambda in searching for optimal lambda
  float_t lambLogMax;
  /// Interval size with logrange used in searching optimal lambda
  float_t lambStepsize;
  /// Number of iteration in searching for optimal lambda
  size_t lambSteps;
  /// Actual base learners number for Adaboosting
  size_t actualBaseLearners;
  /// Judgement of grid refine
  bool refinement;
  /// Select the refine mode(1:use grid number, 2: use grid number percentage)
  size_t refineMode;
  /// Number of refinement with a certain percentage of Grid points
  size_t refineTimes;
  /// Number of Grid points to refine
  size_t numOfAda;
  /// Percentage of Grid points to refine(between 0 and 1)
  float_t perOfAda;
  /// Set the boost mode (1: Discrete Adaboost, 2: Real Adaboost)
  size_t boostMode;

  /**
   * Performs a solver to get alpha
  *
   * @param lambda the regularisation parameter
  * @param weight the weights of examples
  * @param alpha output the coefficients of the sparse grid's basis functions
  * @param final judgement the final step of this base learner
  */
  virtual void alphaSolver(float_t& lambda, SGPP::base::DataVector& weight,
                           SGPP::base::DataVector& alpha, bool final) = 0;

 public:

  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param gridType reference to the of grid type(1 = Linear Grid, 2 = LinearL0Boundary Grid, 3 = ModLinear Grid)
   * @param gridLevel reference to the level of grid
   * @param trainData reference to the training dataset
   * @param trainDataClass reference to the class(real value in regression) of training dataset
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
  AlgorithmAdaBoostBase(SGPP::base::Grid& SparseGrid, size_t gridType,
                        SGPP::base::HashGenerator::level_t gridLevel, SGPP::base::DataMatrix& trainData,
                        SGPP::base::DataVector& trainDataClass, size_t NUM, float_t lambda, size_t IMAX,
                        float_t eps, size_t IMAX_final, float_t eps_final, float_t firstLabel,
                        float_t secondLabel, float_t threshold, float_t maxLambda, float_t minLambda,
                        size_t searchNum, bool refine, size_t refineMode, size_t refineNum,
                        size_t numberOfAda, float_t percentOfAda, size_t mode);


  /**
   * Std-Deconstructor
   */
  virtual ~AlgorithmAdaBoostBase();

  /**
   * Performs the Discrete Adaboost
  *
  * @param hypoWeight the vector to store hypothesis weights(Alpha-t)
  * @param weightError the vector to store the weight error of each iteration
  * @param weights the matrix to store weights of every training date for every weak learner
  * @param decision the matrix to store the decision right or not according to the true class
  * @param testData reference to the testing dataset
  * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
  * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
   */
  void doDiscreteAdaBoost(SGPP::base::DataVector& hypoWeight,
                          SGPP::base::DataVector& weightError, SGPP::base::DataMatrix& weights,
                          SGPP::base::DataMatrix& decision, SGPP::base::DataMatrix& testData,
                          SGPP::base::DataMatrix& algorithmValueTrain,
                          SGPP::base::DataMatrix& algorithmValueTest);

  /**
   * Performs the Real Adaboost
  *
  * @param weights the matrix to store weights of every training date for every weak learner
  * @param testData reference to the testing dataset
  * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
  * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
   */
  void doRealAdaBoost(SGPP::base::DataMatrix& weights,
                      SGPP::base::DataMatrix& testData, SGPP::base::DataMatrix& algorithmValueTrain,
                      SGPP::base::DataMatrix& algorithmValueTest);

  /**
   * Performs the Adaboost.R2(a regression algorithm)
  *
  * @param weights the matrix to store weights of every training date for every weak learner
  * @param testData reference to the testing dataset
  * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
  * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
  * @param lossFucType the loss function type(linear, square or exponential)
   */
  void doAdaBoostR2(SGPP::base::DataMatrix& weights,
                    SGPP::base::DataMatrix& testData, SGPP::base::DataMatrix& algorithmValueTrain,
                    SGPP::base::DataMatrix& algorithmValueTest, std::string lossFucType);

  /**
   * Performs the Adaboost.RT(a regression algorithm)
  *
  * @param weights the matrix to store weights of every training date for every weak learner
  * @param testData reference to the testing dataset
  * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
  * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
  * @param Tvalue the threshold to demarcate the prediction correctness(only from 0 to 1)
  * @param powerType the error rate power coefficient(linear, square or cubic)
   */
  void doAdaBoostRT(SGPP::base::DataMatrix& weights,
                    SGPP::base::DataMatrix& testData, SGPP::base::DataMatrix& algorithmValueTrain,
                    SGPP::base::DataMatrix& algorithmValueTest, float_t Tvalue,
                    std::string powerType);

  /**
   * Performs a real value calculate for the testing dataset
   *
   * @param testData reference to the testing dataset
  * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
  * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
   */
  void eval(SGPP::base::DataMatrix& testData,
            SGPP::base::DataMatrix& algorithmValueTrain,
            SGPP::base::DataMatrix& algorithmValueTest);

  /**
   * Performs a classify for the testing dataset according to the baselearners get from the algorithm
   *
   * @param testData reference to the testing dataset
  * @param algorithmClassTrain reference to the class of training dataset got from the algorithm
  * @param algorithmClassTest reference to the class of testing dataset got from the algorithm
  * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
  * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
   */
  void classif(SGPP::base::DataMatrix& testData,
               SGPP::base::DataVector& algorithmClassTrain,
               SGPP::base::DataVector& algorithmClassTest,
               SGPP::base::DataMatrix& algorithmValueTrain,
               SGPP::base::DataMatrix& algorithmValueTest);

  /**
   * Performs an accuracy evaluation for the testing dataset
   *
   * @param testData reference to the testing dataset
   * @param testDataClass reference to the class of testing dataset
   * @param accuracy_train reference to the accuracy for the training dataset
  * @param accuracy_test reference to the accuracy for the testing dataset
   */
  void getAccuracy(SGPP::base::DataMatrix& testData,
                   SGPP::base::DataVector& testDataClass, float_t* accuracy_train,
                   float_t* accuracy_test);

  /**
   * Performs an evaluation to get ROC related parameter
   *
   * @param validationData reference to the validation dataset
   * @param validationDataClass reference to the class of validation dataset
  * @param acc reference to the accuracy for the validation dataset
  * @param sensitivity reference to the sensitivity for the validation dataset
  * @param specificity reference to the specificity for the validation dataset
  * @param precision reference to the precision for the validation dataset
  * @param recall reference to the recall for the validation dataset
  * @param fOneScore reference to the specificity for the validation dataset
   */
  void getROC(SGPP::base::DataMatrix& validationData,
              SGPP::base::DataVector& validationDataClass, float_t* acc, float_t* sensitivity,
              float_t* specificity, float_t* precision, float_t* recall, float_t* fOneScore);

  /**
   * Performs an accuracy evaluation for the testing dataset with a specified number of base learner
   *
   * @param testData reference to the testing dataset
   * @param testDataClass reference to the class of testing dataset
  * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
  * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
   * @param accuracy_train reference to the accuracy for the training dataset
  * @param accuracy_test reference to the accuracy for the testing dataset
  * @param yourBaseLearner the number of base learner specified
   */
  void getAccuracyBL(SGPP::base::DataMatrix& testData,
                     SGPP::base::DataVector& testDataClass,
                     SGPP::base::DataMatrix& algorithmValueTrain,
                     SGPP::base::DataMatrix& algorithmValueTest, float_t* accuracy_train,
                     float_t* accuracy_test, size_t yourBaseLearner);

  /**
   * Performs refinement of grid to get an adaptive grid
   *
  * @param alpha_ada the coefficients of the sparse grid's basis functions and to be refined
  * @param weight_ada the weights of examples
  * @param curBaseLearner the current base learner
  */
  void doRefinement(SGPP::base::DataVector& alpha_ada,
                    SGPP::base::DataVector& weight_ada, size_t curBaseLearner);

  /**
       * Get the actual base learners after doing adaboosting
       *
  */
  size_t getActualBL();

  /**
       * Get the mean GridPoint ever used in adaboosting
  *
       * @param baseLearner number of baselearner
  */
  size_t getMeanGridPoint(size_t baseLearner);

  /**
       * Get the max GridPoint ever used in adaboosting
       *
       * @param baseLearner number of baselearner
  */
  size_t getMaxGridPoint(size_t baseLearner);

  /**
       * Get the sum GridPoint ever used in adaboosting
       *
       * @param baseLearner number of baselearner
  */
  size_t getSumGridPoint(size_t baseLearner);

  /**
   * Performs a hypothesis classifier
  *
  * @param realValue real value of function
   */
  float_t hValue(float_t realValue);
};
}
}
#endif


