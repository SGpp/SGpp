/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (tjsongzw@gmail.com)
// @author Benjamin (pehersto@in.tum.de)

#include "algorithm/datadriven/AlgorithmAdaBoost.hpp"
#include "exception/operation_exception.hpp"
#include "base/operation/BaseOpFactory.hpp"

namespace sg
{
namespace datadriven
{
	AlgorithmAdaBoost::AlgorithmAdaBoost(sg::base::Grid& SparseGrid, size_t gridType, size_t gridLevel, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, double firstLabel, double secondLabel, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineNum, double percentOfAda, size_t alphaMethod, std::string vecMode)
    {
		if (refine && (gridType != 1 && gridType != 2 && gridType != 3))
		{
			throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid)!");
		}
		if (refine && (percentOfAda >= 1.0 || percentOfAda <= 0.0))
		{
			throw new sg::base::operation_exception("AlgorithmAdaboost : Only number between 0 and 1 is the supported percent to Adaptive!");
		}
		if (alphaMethod != 1 && alphaMethod != 2)
		{
			throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 are supported method to solve alpha(1 : use normal DMWeightMatrix as System Matrix, 2 : use DMWeightMatrixVectorizedIdentity as System Matrix)!");
		}
		if (vecMode != "X86SIMD" && vecMode != "OCL" && vecMode != "ArBB" && vecMode != "HYBRID_X86SIMD_OCL")
		{
			throw new sg::base::operation_exception("AlgorithmAdaboost : Only X86SIMD or OCL or ArBB or HYBRID_X86SIMD_OCL are supported vector extensions!");
		}

		sg::base::GridStorage* gridStorage = SparseGrid.getStorage();
        this->grid = &SparseGrid;
		this->type = gridType;
		this->gridPoint = gridStorage->size();
		maxGridPoint = 0;
		this->level = gridLevel;
        this->lamb = lambda;
        this->data = &trainData;
        this->classes = &trainDataClass;
        this->numData = trainData.getNrows();
        this->dim = gridStorage->dim();
        this->numBaseLearners = NUM;
        this->imax = IMAX;
        this->epsilon = eps;
		this->labelOne = firstLabel;
		this->labelTwo = secondLabel;
		this->lambLogMax = log(maxLambda);
		this->lambSteps = searchNum;
		if(searchNum == 1)
			this->lambStepsize = (log(maxLambda) - log(minLambda))/2;
		else
			this->lambStepsize = (log(maxLambda) - log(minLambda))/(searchNum - 1);
		this->actualBaseLearners = 0;
		this->refinement = refine;
		this->refineTimes = refineNum;
		this->perOfAda = percentOfAda;
		this->alphaMethod = alphaMethod;
		this->vecMode = vecMode;
    }
    
	AlgorithmAdaBoost::~AlgorithmAdaBoost()
	{

	}
	
	void AlgorithmAdaBoost::doAdaBoost(sg::base::DataVector& hypoWeight, sg::base::DataVector& weightError, sg::base::DataMatrix& weights, sg::base::DataMatrix& decision, sg::base::DataMatrix& testData, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest)
	{
		sg::base::DataVector weight(this->numData);
		weight.setAll(1.0/double(this->numData));
		sg::base::OperationEval* opEval = sg::op_factory::createOperationEval(*this->grid);
		// to store certain train data point
		sg::base::DataVector p_train(this->dim);
		// to store certain train data point
		sg::base::DataVector p_test(this->dim);
		// create vector to store the hypothesis of the training data according to certain alpha vector(base learner)
		sg::base::DataVector newclasses(this->numData);
		// create vector to store the identity of comparing hypothesis of training data with class of
		sg::base::DataVector identity(this->numData);
		
		sg::base::DataVector tmpweight(this->numData);

		for (size_t count = 0; count < this->numBaseLearners; count++)
		{
			(this->actualBaseLearners)++;
			std::cout << std::endl;
			std::cout << "This is the " << this->actualBaseLearners << "th weak learner." << std::endl;
			std::cout << std::endl;
			
            // create coefficient vector
			sg::base::DataVector alpha_train(this->gridPoint);
			sg::base::DataVector alpha_learn(this->gridPoint);
			std::cout << "gridPoint: " << this->gridPoint << std::endl;
			if (maxGridPoint < this->gridPoint)
				maxGridPoint = this->gridPoint;
			alpha_train.setAll(0.0);
			weights.setColumn(count, weight);
			
			// calculate alpha
			if (this->alphaMethod == 1)
			{
				sg::base::OperationMatrix* C = sg::op_factory::createOperationIdentity(*this->grid);
				alphaSolver(C, this->lamb, weight, alpha_train);
				delete C;
			}
			else if (this->alphaMethod == 2)
				alphaSolverVectorizedIdentity(this->lamb, weight, alpha_train);
			// should not happen because this exception should have been thrown some lines upwards!
			else
			{
				throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 are supported method to solve alpha(1 : use normal DMWeightMatrix as System Matrix, 2 : use DMWeightMatrixVectorizedIdentity as System Matrix)!");
			}
			   
			if(this->refinement)
			{
				doRefinement(alpha_train, weight);
				opEval = sg::op_factory::createOperationEval(*this->grid);
				alpha_learn.resizeZero(alpha_train.getSize());
			}
			//set the alpha for testing data(copy of alpha for training data)
			for (size_t index = 0; index < alpha_train.getSize(); index++)
				alpha_learn.set(index, alpha_train.get(index));

			for (size_t i = 0; i < this->numData; i++)
			{
				this->data->getRow(i, p_train);
				double value_train = opEval->eval(alpha_train, p_train);
				newclasses.set(i, hValue(value_train));
			}
			
			for (size_t i = 0; i < this->numData; i++)
			{
				if (newclasses.get(i) == this->classes->get(i))
				{
					identity.set(i, 0.0);
					decision.set(i, count, 1.0);
				}
				else
				{
					identity.set(i, 1.0);
					decision.set(i, count, 0.0);
				}
			}
			// calculate the weight error
			weightError.set(count, weight.dotProduct(identity));
			
			// find the optimal lambda to minimize the weighted error
			if (this->lambSteps > 0 && count > 0)
			{
				double cur_lambda;
				double weighterror;
				double minWeightError = weightError.get(count);
				for (size_t it = 0; it < this->lambSteps; it++)
				{
					std::cout << std::endl;
					std::cout << "This is the " << it + 1 << "th search of " << this->actualBaseLearners << "th weak learner." << std::endl;
					std::cout << std::endl;
					alpha_train.setAll(0.0);
					cur_lambda = exp(this->lambLogMax - it*this->lambStepsize);

					if (this->alphaMethod == 1)
					{
						sg::base::OperationMatrix* C_search = sg::op_factory::createOperationIdentity(*this->grid);
						alphaSolver(C_search, cur_lambda, weight, alpha_train);
						delete C_search;
					}
					else if (this->alphaMethod == 2)
						alphaSolverVectorizedIdentity(cur_lambda, weight, alpha_train);
					// should not happen because this exception should have been thrown some lines upwards!
					else
					{
						throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 are supported method to solve alpha(1 : use normal DMWeightMatrix as System Matrix, 2 : use DMWeightMatrixVectorizedIdentity as System Matrix)!");
					}

					for (size_t i = 0; i < this->numData; i++)
					{
						this->data->getRow(i, p_train);
						double value_seach = opEval->eval(alpha_train, p_train);
						newclasses.set(i, hValue(value_seach));
					}
					
					for (size_t i = 0; i < this->numData; i++)
					{
						if (newclasses.get(i) == this->classes->get(i))
						{
							identity.set(i, 0.0);
						}
						else
						{
							identity.set(i, 1.0);
						}
					}

					// compare the weight error we need the minimum weight error 
					weighterror = weight.dotProduct(identity);
					if (weighterror < minWeightError)
					{
						minWeightError = weighterror;
						weightError.set(count, weighterror);
						//reset the alpha for testing data(copy of alpha for training data)
						for (size_t index = 0; index < alpha_train.getSize(); index++)
							alpha_learn.set(index, alpha_train.get(index));
						for (size_t i = 0; i < this->numData; i++)
						{
							if (newclasses.get(i) == this->classes->get(i))
							{
								decision.set(i, count, 1.0);
							}
							else
							{
								decision.set(i, count, 0.0);
							}
						}
					}
				}
			}
			
			// to judge the classification whether match the simple condition of Adaboost
			if (weightError.get(count) >= 0.50)
			{
				std::cout << std::endl << "The training error rate exceeds 0.5 after " << count + 1 << " iterations" << std::endl;
				(this->actualBaseLearners)--;
				break;
			}
			// calculate the weight of this weak classif
			double hypoweight;
			if (weightError.get(count) == 0)
			{
				hypoweight = log(1e+10);
				hypoWeight.set(count, hypoweight);
			}
			else
			{
				hypoweight = 0.5 * (log((1.0 - weightError.get(count))/weightError.get(count)));
				hypoWeight.set(count, hypoweight);
			}

			// calculate the algorithm value of the testing data and training data
			// for training data
			for (size_t i = 0; i < numData; i++)
			{
				this->data->getRow(i, p_train);
				double value_train = opEval->eval(alpha_learn, p_train);
				// when there is only one baselearner actually, we do as following, just use normal classify to get the value
				if (this->numBaseLearners == 1)
					algorithmValueTrain.set(i, count, value_train);
				else if (count == 0)
					algorithmValueTrain.set(i, count, hypoweight * hValue(value_train));
				// each column is the sum value of baselearner respect to the column index
				else
					algorithmValueTrain.set(i, count, algorithmValueTrain.get(i, count - 1) + hypoweight * hValue(value_train));
			}
			// for testing data
			for (size_t i = 0; i < testData.getNrows(); i++)
			{
				testData.getRow(i, p_test);
				double value_test = opEval->eval(alpha_learn, p_test);
				// when there is only one baselearner actually, we do as following, just use normal classify to get the value
				if (this->numBaseLearners == 1)
					algorithmValueTest.set(i, count, value_test);
				else if (count == 0)
					algorithmValueTest.set(i, count, hypoweight * hValue(value_test));
				// each column is the sum value of baselearner respect to the column index
				else
					algorithmValueTest.set(i, count, algorithmValueTest.get(i, count - 1) + hypoweight * hValue(value_test));
			}

			double helper;
			for (size_t i = 0; i < this->numData; i++)
			{
				// helper = weight.get(i)*exp(-hypoWeight.get(count) * newclasses.get(i)*this->classes->get(i));
				if (newclasses.get(i) == this->classes->get(i))
					helper = weight.get(i)*exp(-hypoWeight.get(count));
				else
					helper = weight.get(i)*exp(hypoWeight.get(count));
				tmpweight.set(i, helper);
			}
			// normalization constant, this expression equals to normalizer = 2 * sqrt((weightError.get(count)) * (1.0 - weightError.get(count)));
			double normalizer = tmpweight.sum();
			tmpweight.mult(1.0/normalizer);
			// get new weights vector
			weight = tmpweight;

			if (count < this->numBaseLearners - 1 && this->refinement)
			{
				//reset the grid to the regular grid
				if(this->type == 1)
				{
					this->grid = sg::base::Grid::createLinearGrid(this->dim);
					std::cout << "reset to the regular LinearGrid" << std::endl;
				}
				else if(this->type == 2)
				{
					this->grid = sg::base::Grid::createLinearBoundaryGrid(this->dim);
					std::cout << "reset to the regular LinearBoundaryGrid" << std::endl;
				}
				else if(this->type == 3)
				{
					this->grid = sg::base::Grid::createModLinearGrid(this->dim);
					std::cout << "reset to the regular ModLinearGrid" << std::endl;
				}
				// should not happen because this exception should have been thrown some lines upwards!
				else
				{
					throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 or 3 are supported gridType(1 = Linear Grid, 2 = LinearBoundary Grid, 3 = ModLinear Grid)!");
				}
				sg::base::GridGenerator* gridGen = this->grid->createGridGenerator();
				gridGen->regular(this->level);
				sg::base::GridStorage* gridStorage_reset = this->grid->getStorage();
				std::cout << std::endl;
				std::cout << "# number of grid points:  " << gridStorage_reset->size() << std::endl;
				delete gridGen;
			}
		}
		delete opEval;
	}
	
    void AlgorithmAdaBoost::eval(sg::base::DataMatrix& testData, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest)
    {
		sg::base::DataVector theHypoWeight(this->numBaseLearners);
		sg::base::DataVector theWeightError(this->numBaseLearners);
		sg::base::DataMatrix weightsMatrix(this->numData, this->numBaseLearners);
		sg::base::DataMatrix decisionMatrix(this->numData, this->numBaseLearners);
		theHypoWeight.setAll(0.0);
		theWeightError.setAll(0.0);
		weightsMatrix.setAll(0.0);
		decisionMatrix.setAll(0.0);
		doAdaBoost(theHypoWeight, theWeightError, weightsMatrix, decisionMatrix, testData, algorithmValueTrain, algorithmValueTest);
	}
	
	void AlgorithmAdaBoost::classif(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmClassTrain, sg::base::DataVector& algorithmClassTest, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest)
	{
		eval(testData, algorithmValueTrain, algorithmValueTest);
		for (size_t i = 0; i < this->numData; i++)
		{
			algorithmClassTrain.set(i, hValue(algorithmValueTrain.get(i, this->actualBaseLearners - 1)));
		}
		for (size_t i = 0; i < testData.getNrows(); i++)
		{
			algorithmClassTest.set(i, hValue(algorithmValueTest.get(i, this->actualBaseLearners - 1)));
		}
	}

	void AlgorithmAdaBoost::getAccuracy(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, double* accuracy_train, double* accuracy_test)
	{
		/* get the accuracy */
		size_t right_test = 0;
		size_t right_train = 0;
		sg::base::DataVector classTrain(this->numData);
		sg::base::DataVector classTest(testDataClass.getSize());
		sg::base::DataMatrix valueTrain(this->numData, this->numBaseLearners);
		sg::base::DataMatrix valueTest(testDataClass.getSize(), this->numBaseLearners);
		classif(testData, classTrain, classTest, valueTrain, valueTest);
		// for training data
		for(size_t i = 0; i < this->numData; i++)
		{
			if(classTrain.get(i) == this->classes->get(i))
				right_train = right_train + 1;
			//			std::cout << i << " " << classTrain.get(i) << std::endl;
		}
		*accuracy_train = double(right_train)/double(this->numData);
		// for testing data
		for(size_t i = 0; i < testData.getNrows(); i++)
		{
			if(classTest.get(i) == testDataClass.get(i))
				right_test = right_test + 1;
			//			std::cout << i << " " << classTest.get(i) << std::endl;
		}
		*accuracy_test = double(right_test)/double(classTest.getSize());
	}
	
	void AlgorithmAdaBoost::getAccuracyBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& algorithmValueTrain, sg::base::DataMatrix& algorithmValueTest, double* accuracy_train, double* accuracy_test, size_t yourBaseLearner)
	{
		size_t right_test = 0;
		size_t right_train = 0;

		// for training data
		for (size_t i = 0; i < this->numData; i++)
		{
			if (hValue(algorithmValueTrain.get(i, yourBaseLearner - 1)) == this->classes->get(i))
				right_train = right_train + 1;
		}
		*accuracy_train = double(right_train)/double(this->numData);

		// for testing data
		for (size_t i = 0; i < testData.getNrows(); i++)
		{
			if (hValue(algorithmValueTest.get(i, yourBaseLearner - 1)) == testDataClass.get(i))
				right_test = right_test + 1;
		}
		*accuracy_test = double(right_test)/double(testDataClass.getSize());
	}
	
	void AlgorithmAdaBoost::doRefinement(sg::base::DataVector& alpha_ada, sg::base::DataVector& weight_ada)
	{
		for (size_t adaptiveStep = 1; adaptiveStep <= this->refineTimes; adaptiveStep++)
		{
			sg::base::GridGenerator* myGenerator = this->grid->createGridGenerator();
			size_t refineNumber = this->perOfAda * (this->grid->getSize());
			//force to refine at least one point
			if(refineNumber == 0)
			  refineNumber = 1;

			sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&alpha_ada, refineNumber, 0.0);
			myGenerator->refine(myRefineFunc);
			delete myRefineFunc;
			delete myGenerator;

			sg::base::GridStorage* gridStorage_ada = this->grid->getStorage();
			size_t gridPts = gridStorage_ada->size();
			if (maxGridPoint < gridPts)
				maxGridPoint = gridPts;
			std::cout << std::endl;
			std::cout << "Refinement time step: " << adaptiveStep << ", new grid size: " << gridPts << ", refined number of grid points: " << refineNumber << std::endl;

			// extend alpha vector (new entries uninitialized)
			alpha_ada.resizeZero(gridPts);

			// calculate new alpha

			if (this->alphaMethod == 1)
			{
				sg::base::OperationMatrix* C_ada = sg::op_factory::createOperationIdentity(*this->grid);
				alphaSolver(C_ada, this->lamb, weight_ada, alpha_ada);
				delete C_ada;
			}
			else if (this->alphaMethod == 2)
			{
				alphaSolverVectorizedIdentity(this->lamb, weight_ada, alpha_ada);
			}
			// should not happen because this exception should have been thrown some lines upwards!
			else
			{
				throw new sg::base::operation_exception("AlgorithmAdaboost : Only 1 or 2 are supported method to solve alpha(1 : use normal DMWeightMatrix as System Matrix, 2 : use DMWeightMatrixVectorizedIdentity as System Matrix)!");
			}
		}
	}

	void AlgorithmAdaBoost::alphaSolver(sg::base::OperationMatrix*& C, double& lambda, sg::base::DataVector& weight, sg::base::DataVector& alpha)
	{
		sg::datadriven::DMWeightMatrix WMatrix(*grid, *data, *C, lambda, weight);
		sg::base::DataVector rhs(alpha.getSize());
		WMatrix.generateb(*classes, rhs);
		sg::solver::ConjugateGradients myCG(this->imax, this->epsilon);
		myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
	}

	void AlgorithmAdaBoost::alphaSolverVectorizedIdentity(double& lambda, sg::base::DataVector& weight, sg::base::DataVector& alpha)
	{
		sg::datadriven::DMWeightMatrixVectorizedIdentity WMatrix(*grid, *data, lambda, weight, this->vecMode);
		sg::base::DataVector rhs(alpha.getSize());
		WMatrix.generateb(*classes, rhs);
		sg::solver::ConjugateGradients myCG(this->imax, this->epsilon);
		myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
	}


	size_t AlgorithmAdaBoost::getActualBL()
	{
		return this->actualBaseLearners;
	}

	size_t AlgorithmAdaBoost::getGridPoint()
	{
		return maxGridPoint;
	}
}
}
