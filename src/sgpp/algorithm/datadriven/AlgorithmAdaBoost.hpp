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
#include "algorithm/datadriven/DMWeightMatrixVectorizedIdentity.hpp"
#include "operation/common/OperationHierarchisation.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include <math.h>
#include <string>
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
        size_t dim;
            /// Pointer to the class of the data vector
		sg::base::DataVector* classes;
            /// Number of base learner for Adaboosting
        size_t numBaseLearners;
		    /// the grid
		sg::base::Grid* grid;
		    /// type of grid possible value are 1, 2 or 3(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid);
		size_t type;
    		/// Number of grid points
		size_t gridPoint;
		    /// grid level
		size_t level;
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
		    /// Judgement of grid refine
		bool refinement;
    		/// Number of refinement with a certain percentage of Grid points
		size_t refineTimes;
		    /// Percentage of Grid points to refine(between 0 and 1)
		double perOfAda;
		    /// Different SystemMatrixWithWeights to solve alpha, possible value are 1 or 2(1 = DMWeightMatrix, 2 = DMWeightMatrixVectorizedIdentity)
		size_t alphaMethod;
    		/// Vectorization mode, possible values are SSE, AVX, OCL, ArBB
		std::string vecMode;


    public:

            /**
             * Std-Constructor
             * 
             * @param SparseGrid reference to the sparse grid
             * @param gridType reference to the of grid type(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid)
             * @param gridLevel reference to the level of grid
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
			 * @param refine the judgement of refine
			 * @param refineNum the Number of refinement with a certain percentage of Grid points
			 * @param percentOfAda the percentage of Grid points to refine
			 * @param alphaMethod different SystemMatrixWithWeights to solve alpha, possible value are 1 or 2(1 = DMWeightMatrix, 2 = DMWeightMatrixVectorizedIdentity)
			 * @param vecMode vectorization mode, possible values are SSE, AVX, OCL, ArBB
             */
        AlgorithmAdaBoost(sg::base::Grid& SparseGrid, size_t gridType, size_t gridLevel, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, double firstLabel, double secondLabel, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineNum, double percentOfAda, size_t alphaMethod, std::string vecMode);
        

            /**
             * Std-Deconstructor
             */
        virtual ~AlgorithmAdaBoost();

            /**
             * Performs the algorithm
			 *
			 * @param hypoWeight the vector to store hypothesis weights(Alpha-t)
			 * @param weightError the vector to store the weight error of each iteration
			 * @param weights the matrix to store weights of every training date for every weak learner
			 * @param decision the matrix to store the decision right or not according to the true class
			 * @param testData reference to the testing dataset
			 * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
			 * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
             */
        void doAdaBoost(sg::base::DataVector& hypoWeight, sg::base::DataVector& weightError, sg::base::DataMatrix& weights, sg::base::DataMatrix& decision, sg::base::DataMatrix& testData, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest);

            /**
             * Performs a real value calculate for the testing dataset
             * 
             * @param testData reference to the testing dataset
			 * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
			 * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
             */
        void eval(sg::base::DataMatrix& testData,  sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest);

            /**
             * Performs a classify for the testing dataset according to the baselearners get from the algorithm
             * 
             * @param testData reference to the testing dataset
			 * @param algorithmClassTrain reference to the class of training dataset got from the algorithm
			 * @param algorithmClassTest reference to the class of testing dataset got from the algorithm
			 * @param algorithmValueTrain the matrix reference to the real value of training dataset got from the algorithm with diff base learners
			 * @param algorithmValueTest the matrix reference to the real value of testing dataset got from the algorithm with diff base learners
             */
        void classif(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmClassTrain, sg::base::DataVector& algorithmClassTest, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest);
        
            /**
             * Performs an accuracy evaluation for the testing dataset
             *
             * @param testData reference to the testing dataset
             * @param testDataClass reference to the class of testing dataset
             * @param accuracy_train reference to the accuracy for the training dataset
			 * @param accuracy_test reference to the accuracy for the testing dataset
             */
        void getAccuracy(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, double* accuracy_train, double* accuracy_test);

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
		void getAccuracyBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest, double* accuracy_train, double* accuracy_test, size_t yourBaseLearner);

            /**
             * Performs refinement of grid to get an adaptive grid
             *
			 * @param alpha_ada the coefficients of the sparse grid's basis functions and to be refined
			 * @param weight_ada the weights of examples
			 */
		void doRefinement(sg::base::DataVector& alpha_ada, sg::base::DataVector& weight_ada);

            /**
             * Performs a solver to get alpha use DMWeightMatrix as the System Matrix
			 *
             * @param C OperationMatrix for the regularisation mehtod
             * @param lambda the regularisation parameter
			 * @param weight the weights of examples
			 * @param alpha output the coefficients of the sparse grid's basis functions
			 */
		void alphaSolver(sg::base::OperationMatrix*& C, double& lambda, sg::base::DataVector& weight, sg::base::DataVector& alpha);

            /**
             * Performs a solver to get alpha use DMWeightMatrixVectorizedIdentity as the System Matrix
			 *
             * @param lambda the regularisation parameter
			 * @param weight the weights of examples
			 * @param alpha output the coefficients of the sparse grid's basis functions
			 */
		void alphaSolverVectorizedIdentity(double& lambda, sg::base::DataVector& weight, sg::base::DataVector& alpha);

		    /**
             * Get the actual base learners after doing adaboosting
             *
			 */
		size_t getActualBL();
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
		
	};
}
}
#endif

        
        
