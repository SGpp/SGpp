/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (tjsongzw@gmail.com)
// @author Benjamin (pehersto@in.tum.de)
#ifndef ALGORITHMADABOOST_HPP
#define ALGORITHMADABOOST_HPP

#include "grid/GridStorage.hpp"
#include "grid/Grid.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "operation/datadriven/OperationMultipleEval.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "algorithm/datadriven/DMWeightMatrix.hpp"
#include <math.h>
#include <vector>
#include <utility>
#include <iostream>
#include <cstdlib>

namespace sg
{
namespace datadriven
{

/*
 * Algorithm for Adaboosting
 * This algorithm is to train Train base learner according to sample distribution and obtain hypothesis
 * get the hypothesis weight
 * Then combine hypothesis linearly
 *
 * The main idea behind the algorithm is to get a better accuracy in classify dateset according to the training dataset
 *
 */
    class AlgorithmAdaBoost
    {

    private:
            /// the lambda, the regularisation parameter
        double lamb;
            /// the size of the grid
        size_t numData;
            /// Pointer to the data matrix
		sg::base::DataMatrix* data;
            /// DataSet Dimension
        int dim;
            /// Pointer to the class of the data vector
		sg::base::DataVector* classes;
            /// Number of base learner for Adaboosting
        size_t numBaseLearners;
		    /// type of the grid
		sg::base::Grid* grid;
        	/// OperationMatrix, the regularisation mehtod
		sg::base::OperationMatrix* C;
            /// Parameter for CG solver
        size_t imax;
		    /// Parameter for CG solver
        double epsilon;
    		/// One label of the DataSet
		double labelOne;
		    /// Another label of the DataSet
		double labelTwo;
    		/// Log of the Max lambda in searching for optimal lambda
		double lambLogMax;
     		/// Interval size with logrange used in searching optimal lambda
		double lambStepsize;
    		/// Number of iteration in searching for optimal lambda 
		size_t lambSteps;
    		/// Actual base learners number for Adaboosting
		size_t actualBaseLearners;
            /**
             * Performs a hypothesis classifier
			 *
			 * @param realvalue real value of function 
             */
		double hValue(double realValue)
		{
			double meanValue = (this->labelOne + this->labelTwo) / 2;
			if (realValue > meanValue)
			{
				if (labelOne > labelTwo)
					return labelOne;
				else 
					return labelTwo;
			}
			else 
			{
				if (labelOne > labelTwo)
					return labelTwo;
				else 
					return labelOne;
			}
		};

    public:

            /**
             * Std-Constructor
             * 
             * @param SparseGrid reference to the sparse grid
             * @param trainData reference to the training dataset
             * @param trainDataClass reference to the class of training dataset
             * @param NUM the number of baselearner for Adaboosting
             * @param lambda the regularisation parameter
             * @param IMAX the parameter for ConjugateGradients
             * @param eps the parameter for ConjugateGradients
             * @param firstLabel one label from training dataset
			 * @param secondLabel another label from training dataset
			 * @param maxLambda the max lambda used in searching optimal lambda
			 * @param minLambda the min lambda used in searching optimal lambda
			 * @param searchNum the searching times used in searching for optimal lambda
             */
        AlgorithmAdaBoost(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, double firstLabel, double secondLabel, double maxLambda, double minLambda, size_t searchNum);
        

            /**
             * Std-Deconstructor
             */
        virtual ~AlgorithmAdaBoost();

            /**
             * Performs the algorithm
			 *
			 * @param storageAlpha the matrix to store alpha for each different weights
			 * @param hypoWeight the vector to store hypothesis weights(Alpha-t)
			 * @param weightError the vector to store the weight error of each iteration
			 * @param weights the matrix to store weights of every training date for every weak learner
			 * @param decision the matrix to store the decision right or not according to the true class
             */
        void doAdaBoost(sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, sg::base::DataVector& weightError, sg::base::DataMatrix& weights, sg::base::DataMatrix& decision);

            /**
             * Performs a real value calculate for the testing dataset
             * 
             * @param testData reference to the testing dataset
             * @param algorithmValue reference to the real value got from the algorithm
             */
        void eval(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmValue);

            /**
             * Performs a classify for the testing dataset according to the baselearners get from the algorithm
             * 
             * @param testData reference to the testing dataset
			 * @param algorithmClass reference to the class got from the algorithm
			 * @param algorithmValue reference to the function real value got from the algorithm 
             */
        void classif(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmClass, sg::base::DataVector& algorithmValue);
        
            /**
             * Performs an accuracy and interpolation error evaluation for the testing dataset
			 *
             * @param testData reference to the testing dataset
             * @param testDataClass reference to the class of testing dataset
			 * @param accruracy reference to the the accruracy 
			 * @param error reference to the interpolation error
             */
        void getAccuracyAndError(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, double* accruracy, double* error);

            /**
             * Performs an accuracy evaluation for the testing dataset
             *
             * @param testData reference to the testing dataset
             * @param testDataClass reference to the class of testing dataset
             */
        double getAccuracy(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass);

            /**
             * Performs an interpolation error evaluation for the testing dataset
             *
             * @param testData reference to the testing dataset   
             * @param testDataClass reference to the class of testing dataset
             */
        double getError(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass);

            /**
             * Performs a real value calculate for the testing dataset with a specified number of base learner
             * 
             * @param testData reference to the testing dataset
             * @param algorithmValue reference to the real value got from the algorithm
			 * @param storageAlpha the matrix to store alpha for each different weights
			 * @param hypoWeight the vector to store hypothesis weights(Alpha-t)
			 * @param yourBaseLearner the number of base learner specified
             */
		void evalBL(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmValue, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t yourBaseLearner);

            /**
             * Performs a real value calculate for the testing dataset with a specified number of base learner
             * 
             * @param testData reference to the testing dataset
             * @param algorithmClass reference to the class got from the algorithm
             * @param algorithmValue reference to the real value got from the algorithm
			 * @param storageAlpha the matrix to store alpha for each different weights
			 * @param hypoWeight the vector to store hypothesis weights(Alpha-t)
			 * @param yourBaseLearner the number of base learner specified
             */
		void classifBL(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmClass, sg::base::DataVector& algorithmValue, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t yourBaseLearner);

            /**
             * Performs an accuracy and interpolation error evaluation for the testing dataset with a specified number of base learner
			 *
             * @param testData reference to the testing dataset
             * @param testDataClass reference to the class of testing dataset
			 * @param storageAlpha the matrix to store alpha for each different weights
			 * @param hypoWeight the vector to store hypothesis weights(Alpha-t)
			 * @param accruracy reference to the the accruracy 
			 * @param error reference to the interpolation error
			 * @param yourBaseLearner the number of base learner specified
             */
		void getAccuracyAndErrorBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, double* accruracy, double* error, size_t yourBaseLearner);

            /**
             * Performs an accuracy evaluation for the testing dataset with a specified number of base learner
             *
             * @param testData reference to the testing dataset
             * @param testDataClass reference to the class of testing dataset
			 * @param storageAlpha the matrix to store alpha for each different weights
			 * @param hypoWeight the vector to store hypothesis weights(Alpha-t)
			 * @param yourBaseLearner the number of base learner specified
             */
		double getAccuracyBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t yourBaseLearner);

            /**
             * Performs an interpolation error evaluation for the testing dataset
             *
             * @param testData reference to the testing dataset
             * @param testDataClass reference to the class of testing dataset
			 * @param storageAlpha the matrix to store alpha for each different weights
			 * @param hypoWeight the vector to store hypothesis weights(Alpha-t)
			 * @param yourBaseLearner the number of base learner specified
             */
		double getErrorBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t yourBaseLearner);
    
		    /**
             * Get the actual base learners after doing adaboosting
             *
			 */
		size_t getActualBL();
	};
    
}
}
#endif

        
        
