// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSGD_HPP
#define LEARNERSGD_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
//#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

#include <sgpp/globaldef.hpp>

#include <list>
#include <string>
#include <vector>


namespace sgpp {
namespace datadriven {

class LearnerSGD {
 public:
  LearnerSGD(sgpp::base::RegularGridConfiguration& gridConfig,
             sgpp::base::AdpativityConfiguration& adaptivityConfig
             /*sgpp::datadriven::RegularizationType& regularization*/);
  
  /*
   * Initializes SGD learner.
   *
   * @param pTrainData training dataset: x values
   * @param pTrainLabels training dataset: y values
   * @param pTestData
   * @param pTestLabels
   * @param lambda regularization factor
   * @param gamma step width
   * @param smoothedErrorDecline
   * @param batchSize
   * @param bufferSize
   * */
  virtual void initialize(sgpp::base::DataMatrix& pTrainData,
                          sgpp::base::DataVector& pTrainLabels,
                          sgpp::base::DataMatrix& pTestData,
                          sgpp::base::DataVector& pTestLabels,
                          double lambda,
                          double gamma,
                          double smoothedErrorDecline,
                          size_t batchSize,
                          size_t bufferSize);

  /*
   * Implements stochastic gradient descent.
   * */
  virtual void train(size_t dataNum);

  virtual double getAccuracy(sgpp::base::DataMatrix& testData, sgpp::base::DataVector& testLabels, 
                             double threshold);

  virtual void storeResults(base::DataMatrix& testDataset,
                            base::DataVector& testLabels,
                            double threshold);

  virtual ~LearnerSGD();

  double error;

 protected:

 private:
  /**
   * Generates a regular grid.
   * @return grid
   */
  std::shared_ptr<base::Grid> createRegularGrid();

  /**
   * Computes speciefied error type (e.g. MSE).
   * @param
   * @param
   * @param
   * @return 
   */
  double getError(sgpp::base::DataMatrix& data, sgpp::base::DataVector& labels, 
                  std::string errorType);

  /**
   * Computes error .
   * @param
   * @param
   * @param
   */
  void getBatchError(sgpp::base::DataMatrix& data, sgpp::base::DataVector& labels, 
                     sgpp::base::DataVector& error);

  double getAccuracy(sgpp::base::DataVector& computedLabels,
                     sgpp::base::DataVector& testLabels,
                     double threshold);

  void predict(sgpp::base::DataMatrix& testData,
               sgpp::base::DataVector& computedLabels);

  /**
   * Stores currently processed data point.
   * @param
   * @param
   */
  void pushToBatch(sgpp::base::DataVector& x, double y);
  
  //int getRandom(int limit);

  std::shared_ptr<base::Grid> grid;
  std::shared_ptr<base::DataVector> alpha;
  std::shared_ptr<base::DataVector> alphaAvg;
  std::shared_ptr<base::DataMatrix> trainData;
  std::shared_ptr<base::DataVector> trainLabels;
  std::shared_ptr<base::DataMatrix> testData;
  std::shared_ptr<base::DataVector> testLabels;
  std::shared_ptr<base::DataMatrix> batchData;
  std::shared_ptr<base::DataVector> batchLabels;
  std::shared_ptr<base::DataVector> batchError;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::base::AdpativityConfiguration adaptivityConfig;
  //sgpp::datadriven::RegularizationType regularization;

  std::list<double> smoothedErrorDeclineBuffer;
  double currentBatchError;
  double lambda;
  double gamma;
  double currentGamma;
  double smoothedErrorDecline;

  size_t batchSize;
  size_t smoothedErrorDeclineBufferSize;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERSGD_HPP */

