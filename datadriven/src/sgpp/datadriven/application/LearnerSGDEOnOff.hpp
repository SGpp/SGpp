// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERSGDEONOFF_HPP
#define LEARNERSGDEONOFF_HPP

//#include "Classification.hpp"
#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
//#include "density/Density.hpp"
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

using namespace std;

//namespace sgpp {
//namespace datadriven {

class LearnerSGDEOnOff : public DBMatOnline {
public:
  LearnerSGDEOnOff(sgpp::datadriven::DBMatDensityConfiguration& dconf,
                   sgpp::base::DataMatrix& trainData, sgpp::base::DataVector& trainData_c,
                   sgpp::base::DataMatrix& testData, sgpp::base::DataVector& testData_c,
                   /*double* classLabels, */int classNumber,
                   double lambda, bool usePrior = true, double beta = 0.);
  /*virtual */~LearnerSGDEOnOff();

  /**
   * Trains the learner with the given dataset
   *
   * @param batch_size size of subset of samples used for each training step
   * @param next_cv_step determines when next cross validation has to be triggered
   */
  virtual void train(size_t batch_size, unsigned int next_cv_step, size_t dataNum);
	
	/**
	 * Trains the learner with the given dataset
	 *
	 * @param data the training dataset
	 * @param trainData classes corresponding to the training dataset
	 * @param usePrior use prior density to compute posterior 
	 * @param RefineCoarse_ vector of a pair of a list representing indices of removed grid points and an unsigned int representing added grid points. The vector is of length 'number of classes'
	 */
	virtual void train(sgpp::base::DataMatrix& trainData, sgpp::base::DataVector& classes, bool do_cv = false, std::vector<std::pair<std::list<size_t>, unsigned int> >*  RefineCoarse_ = NULL);
	
	/**
	 * Trains the learner with a dataset that is already split up into its different classes
	 *
	 * @param trainDataClasses A list of pairs. Each pair contains the data points that belong to one class and the corresponding class label
	 * @param use prior density to compute posteriori
	 * @param RefineCoarse_ vector of a pair of a list representing indices of removed grid points and an unsigned int representing added grid points. The vector is of length 'number of classes'
	 */
	virtual void train(std::vector<std::pair<sgpp::base::DataMatrix*, double> >& trainDataClasses, bool do_cv = false, std::vector<std::pair<std::list<size_t>, unsigned int> >*  RefineCoarse_ = NULL);

        virtual double getAccuracy();
	
	virtual sgpp::base::DataVector predict(sgpp::base::DataMatrix* test);
	
	virtual int predict(sgpp::base::DataVector &p);

        virtual void storeResults();
	
	/**
	 * Returns the values of all density functions for one point
	 * Can be used to plot the density functions or to compute confidence values
	 *
	 * @param point The point for which the density functions should be evaluated
	 */
	virtual sgpp::base::DataVector getDensities(sgpp::base::DataVector& point);
	
	/**
	 * Sets the crossvalidation parameters.
	 * They get directly passed onto the DBMatOnlineDE class-instance
	 * 
	 * @param lambda_step how many different lambdas are tried out
	 * @param lambda_start The smallest possible lambda
	 * @param lambda_end The biggest possible lambda
	 * @param test The test matrix
	 * @param test_cc The results of the points in the test matrix
	 * @param logscale Indicates whether the values between lambda_start and lambda_end are searched using logscale or not
	 */
	virtual void setCrossValidationParameters(int lambda_step, double lambda_start, double lambda_end, sgpp::base::DataMatrix *test, sgpp::base::DataMatrix *test_cc, bool logscale);
	
	/**
	 * In case of crossvalidation, returns the current best lambda
	 */
	virtual double getBestLambda();
	
	//virtual void init(std::map<double, int> entriesPerClass);
        virtual void init();
	
	virtual unsigned int getNumClasses();
	
	//void normalize(size_t samples = 1000);
	
	//Get density functions mapped to class labels
	virtual std::vector<std::pair<DBMatOnlineDE*, double> >* getDestFunctions();

	time_t traintime;
	std::map<double, double> prior;

        double error;

protected:  
        base::DataMatrix* trainData;
	base::DataVector* trainLabels;
        //std::shared_ptr<base::DataMatrix> trainData;
        //std::shared_ptr<base::DataVector> trainLabels;
        base::DataMatrix* testData;
        base::DataVector* testLabels;

        DBMatOffline* offline;
	double* classLabels;
	int classNumber;
	std::vector<std::pair<DBMatOnlineDE*, double> >* destFunctions_;
	bool trained_;
	bool usePrior_;
	size_t processedPoints;
	bool initDone;
	double beta_;

	int cv_save_lambda_step;
	double cv_save_lambda_start, cv_save_lambda_end;
	bool cv_save_logscale, cv_saved;
	sgpp::base::DataMatrix *cv_save_test, *cv_save_test_cc;
};

//}  // namespace datadriven
//}  // namespace sgpp

#endif /* LEARNERSGDEONOFF_HPP */
