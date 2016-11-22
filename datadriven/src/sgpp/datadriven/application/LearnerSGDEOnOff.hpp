// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSGDEONOFF_HPP
#define LEARNERSGDEONOFF_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

namespace sgpp {
namespace datadriven {

class LearnerSGDEOnOff : public DBMatOnline {
 public:

  LearnerSGDEOnOff(sgpp::datadriven::DBMatDensityConfiguration& dconf,
                   base::DataMatrix& trainData, base::DataVector& trainDataLabels,
                   base::DataMatrix& testData, base::DataVector& testDataLabels,
                   base::DataMatrix* validData, base::DataVector* validDataLabels,
                   base::DataVector& classLabels, int classNumber,
                   bool usePrior, double beta, double lambda);

  virtual ~LearnerSGDEOnOff();

  /**
   * Trains the learner with the given dataset
   *
   * @param batchSize size of subset of samples used for each training step
   * @param nextCvStep determines when next cross validation has to be triggered
   */
  void train(size_t batchSize, size_t maxDataPasses, 
             string refType, string refMonitor, 
             size_t refPeriod, double accDeclineThreshold, 
             size_t accDeclineBufferSize, size_t minRefInterval,
             bool enableCv, unsigned int nextCvStep, size_t dataNum);
	
  /**
   * Trains the learner with the given dataset
   *
   * @param data the training dataset
   * @param trainData classes corresponding to the training dataset
   * @param usePrior use prior density to compute posterior 
   * @param RefineCoarse_ vector of a pair of a list representing indices 
   *        of removed grid points and an unsigned int representing added grid points.
   */
  void train(base::DataMatrix& trainData, base::DataVector& classes, bool do_cv = false, 
             std::vector<std::pair<std::list<size_t>, unsigned int> >* RefineCoarse_ = nullptr);
	
  /**
   * Trains the learner with a dataset that is already split up into its different classes
   *
   * @param trainDataClasses A list of pairs. Each pair contains the data points that belong to 
   *        one class and the corresponding class label.
   * @param use prior density to compute posteriori
   * @param RefineCoarse_ vector of a pair of a list representing indices of 
   *        removed grid points and an unsigned int representing added grid points.
   */
  void train(std::vector<std::pair<base::DataMatrix*, double> >& trainDataClasses, bool do_cv = false, 
             std::vector<std::pair<std::list<size_t>, unsigned int> >*  RefineCoarse_ = nullptr);

  double getAccuracy();
	
  sgpp::base::DataVector predict(base::DataMatrix* test);
	
  int predict(base::DataVector &p);

  double getError(base::DataMatrix* data, base::DataVector* labels, 
                  std::string errorType);

  void storeResults();
	
  /**
   * Returns the values of all density functions for one point
   * Can be used to plot the density functions or to compute confidence values
   *
   * @param point The point for which the density functions should be evaluated
   */
  base::DataVector getDensities(base::DataVector& point);
	
  /**
   * Sets the crossvalidation parameters.
   * They get directly passed onto the DBMatOnlineDE class-instance
   * 
   * @param lambda_step how many different lambdas are tried out
   * @param lambda_start The smallest possible lambda
   * @param lambda_end The biggest possible lambda
   * @param test The test matrix
   * @param test_cc The results of the points in the test matrix
   * @param logscale Indicates whether the values between lambda_start 
   *        and lambda_end are searched using logscale or not
   */
  void setCrossValidationParameters(int lambda_step, double lambda_start, 
                                    double lambda_end, base::DataMatrix *test, 
                                    base::DataMatrix *test_cc, bool logscale);
	
  /**
  * In case of crossvalidation, returns the current best lambda
  */
  double getBestLambda();
	
  void init();
	
  unsigned int getNumClasses();
	
  //Get density functions mapped to class labels
  std::vector<std::pair<DBMatOnlineDE*, double> >* getDestFunctions();

  //time_t traintime;
  std::map<double, double> prior;

  double error;
  base::DataVector avgErrors;

 protected:  
  base::DataMatrix* trainData;
  base::DataVector* trainLabels;
  base::DataMatrix* testData;
  base::DataVector* testLabels;
  base::DataMatrix* validData;
  base::DataVector* validLabels;

  base::DataVector classLabels;
  int classNumber;
  std::vector<std::pair<DBMatOnlineDE*, double> >* destFunctions_;
  bool trained;
  bool initDone;
  bool usePrior;	
  double beta;

  DBMatOffline* offline;
  size_t processedPoints;

  int cv_save_lambda_step;
  double cv_save_lambda_start, cv_save_lambda_end;
  bool cv_save_logscale, cv_saved;
  base::DataMatrix *cv_save_test, *cv_save_test_cc;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* LEARNERSGDEONOFF_HPP */
