// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSVM_HPP
#define LEARNERSVM_HPP

#include <sgpp/datadriven/application/PrimalDualSVM.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace datadriven {

class LearnerSVM {
 protected:
  std::shared_ptr<base::Grid> grid;
  std::shared_ptr<base::DataMatrix> trainData;
  std::shared_ptr<base::DataVector> trainLabels;
  std::shared_ptr<base::DataMatrix> testData;
  std::shared_ptr<base::DataVector> testLabels;
  std::shared_ptr<base::DataMatrix> validData;
  std::shared_ptr<base::DataVector> validLabels;
  
  std::shared_ptr<PrimalDualSVM> svm;

  sgpp::base::RegularGridConfiguration gridConfig;
  //sgpp::base::AdpativityConfiguration adaptivityConfig;
  
  /**
   * generates a regular grid
   * @return grid
   */
  std::shared_ptr<base::Grid> createRegularGrid();
  
 public:
   /**
   * Constructor
   *
   * @param gridConfig grid configuration
   */
  LearnerSVM(sgpp::base::RegularGridConfiguration& gridConfig);
  
  virtual ~LearnerSVM();
  
  virtual void initialize(base::DataMatrix& pTrainData, base::DataVector& pTrainLabels,
                          base::DataMatrix& pTestData, base::DataVector& pTestLabels,
                          std::shared_ptr<base::DataMatrix> pValidData,
                          std::shared_ptr<base::DataVector> pValidLabels);

  /**
   * Implements support vector learning with sparse grid kernels.
   *
   * @param 
   * */
  virtual void train(size_t dataNum);

  virtual void storeResults(sgpp::base::DataMatrix& testDataset,
                            sgpp::base::DataVector& testLabels,
                            double threshold);

  virtual double getAccuracy(sgpp::base::DataMatrix& testDataset,
                             const sgpp::base::DataVector& classesReference,
                             const double threshold);

  virtual double getAccuracy(const sgpp::base::DataVector& classesComputed,
                             const sgpp::base::DataVector& classesReference,
                             const double threshold);

  virtual void predict(sgpp::base::DataMatrix& testData,
                       sgpp::base::DataVector& computedLabels);  

  virtual double getError(sgpp::base::DataMatrix& data, sgpp::base::DataVector& labels, 
                          std::string errorType);

  double error;
  sgpp::base::DataVector avgErrors;

};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERSVM_HPP */

