/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (tjsongzw@gmail.com)
// @author Benjamin (pehersto@in.tum.de)

#include "algorithm/datadriven/AlgorithmAdaBoost.hpp"
#include "exception/algorithm_exception.hpp"
#include "base/operation/BaseOpFactory.hpp"

namespace sg
{
namespace datadriven
{
	AlgorithmAdaBoost::AlgorithmAdaBoost(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, double firstLabel, double secondLabel, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineNum)
    {
		sg::base::GridStorage* gridStorage = SparseGrid.getStorage();
        this->grid = &SparseGrid;
		this->gridPoint = gridStorage->size();
        this->lamb = lambda;
        this->data = &trainData;
        this->classes = &trainDataClass;
        this->numData = trainData.getNrows();
        this->dim = gridStorage->dim();
        this->numBaseLearners = NUM;
		this->C = sg::op_factory::createOperationIdentity(SparseGrid);
        this->imax = IMAX;
        this->epsilon = eps;
		this->labelOne = firstLabel;
		this->labelTwo = secondLabel;
		this->lambLogMax = log(maxLambda); 
		this->lambStepsize = (log(maxLambda) - log(minLambda))/(searchNum - 1);
		this->lambSteps = searchNum;
		this->actualBaseLearners = 0;
		this->refinement = refine;
		this->refineTimes = refineNum;
    }
    
	AlgorithmAdaBoost::~AlgorithmAdaBoost()
	{
		delete this->C;
	}
	
	void AlgorithmAdaBoost::doAdaBoost(sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, sg::base::DataVector& weightError, sg::base::DataMatrix& weights, sg::base::DataMatrix& decision)
	{
		sg::base::DataVector weight(this->numData);
		weight.setAll(1.0/double(this->numData));
		sg::base::OperationEval* opEval = sg::op_factory::createOperationEval(*this->grid);
		// to store certain train data point
		sg::base::DataVector p(this->dim);
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
			sg::base::DataVector alpha(this->gridPoint);
			sg::base::DataVector rhs(this->gridPoint);
			alpha.setAll(0.0);
			rhs.setAll(0.0);
			weights.setColumn(count, weight);
			
			sg::datadriven::DMWeightMatrix WMatrix(*grid, *data, *C, this->lamb, weight);
			WMatrix.generateb(*classes, rhs);
			sg::solver::ConjugateGradients myCG(this->imax, this->epsilon);
			myCG.solve(WMatrix, alpha, rhs, false, true, -1.0);

			if(count == 0 && this->refinement)
			{
				for (size_t adaptiveStep = 1; adaptiveStep <= this->refineTimes; adaptiveStep++)
					doRefinement(alpha, rhs, weight, adaptiveStep);
				opEval = sg::op_factory::createOperationEval(*this->grid);
				// resize the DataMatrix storageAlpha
				storageAlpha.resizeZero(alpha.getSize());
			}

			storageAlpha.setColumn(count, alpha);
			
			for (size_t i = 0; i < this->numData; i++)
			{
				this->data->getRow(i, p);
				double value01 = opEval->eval(alpha, p);
				newclasses.set(i, hValue(value01));
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
			if (count > 0)
			{
				double cur_lambda;
				double weighterror;
				double minWeightError = weightError.get(count);
				for (size_t it = 0; it < this->lambSteps; it++)
				{
					std::cout << std::endl;
					std::cout << "This is the " << it + 1 << "th search of " << this->actualBaseLearners << "th weak learner." << std::endl;
					std::cout << std::endl;
					alpha.setAll(0.0);
					cur_lambda = exp(this->lambLogMax - it*this->lambStepsize);
					
					sg::datadriven::DMWeightMatrix WMatrix(*grid, *data, *C, cur_lambda, weight);
					WMatrix.generateb(*classes, rhs);
					sg::solver::ConjugateGradients myCG(this->imax, this->epsilon);
					myCG.solve(WMatrix, alpha, rhs, false, true, -1.0);
					
					for (size_t i = 0; i < this->numData; i++)
					{
						this->data->getRow(i, p);
						double value11 = opEval->eval(alpha, p);
						newclasses.set(i, hValue(value11));
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
						storageAlpha.setColumn(count, alpha);
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
			if (weightError.get(count) == 0)
			{
				std::cout << std::endl << "The training error rate is 0 after  " << count + 1 << " iterations" << std::endl;
				(this->actualBaseLearners)--;
				break;
			}
			
			// calculate the weight of this weak classif
			hypoWeight.set(count, 0.5 * (log((1.0 - weightError.get(count))/weightError.get(count))));
			
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
		}
	}
	
    void AlgorithmAdaBoost::eval(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmValue)
    {
		sg::base::DataMatrix storageOfAlpha(this->gridPoint, this->numBaseLearners);
		sg::base::DataVector theHypoWeight(this->numBaseLearners);
		sg::base::DataVector theWeightError(this->numBaseLearners);
		sg::base::DataMatrix weightsMatrix(this->numData, this->numBaseLearners);
		sg::base::DataMatrix decisionMatrix(this->numData, this->numBaseLearners);
		storageOfAlpha.setAll(0.0);
		theHypoWeight.setAll(0.0);
		theWeightError.setAll(0.0);
		weightsMatrix.setAll(0.0);
		decisionMatrix.setAll(0.0);
		
		doAdaBoost(storageOfAlpha, theHypoWeight, theWeightError, weightsMatrix, decisionMatrix);
		
		sg::base::DataMatrix baseLearnerMatr(testData.getNrows(), this->actualBaseLearners);
		sg::base::DataVector p(this->dim);
		sg::base::OperationEval* operateEval = sg::op_factory::createOperationEval(*this->grid);
		sg::base::DataVector alphaT(storageOfAlpha.getNrows());
		sg::base::DataVector baseLearnerVec(testData.getNrows());
		for (size_t t = 0; t < this->actualBaseLearners; t++)
		{
			storageOfAlpha.getColumn(t, alphaT);
			baseLearnerVec.setAll(0.0);
			for (size_t i = 0; i < testData.getNrows(); i++)
			{
				testData.getRow(i, p);
				double value01 = operateEval->eval(alphaT, p);
				// when there is only one baselearner actually, we do as following, just use normal classify to get the value
				if (this->actualBaseLearners == 1)
					algorithmValue.set(i, value01);
				else
					baseLearnerVec.set(i, hValue(value01));
			}
			if (this->actualBaseLearners > 1)
				baseLearnerMatr.setColumn(t, baseLearnerVec);
		}
		if (this->actualBaseLearners > 1)
		{
			double value02;
			sg::base::DataVector tmp(this->actualBaseLearners);
			for (size_t i = 0; i < testData.getNrows(); i++)
			{
				value02 = 0.0;
				for(size_t k = 0; k < this->actualBaseLearners; k++)
					tmp.set(k, baseLearnerMatr.get(i, k));
				for(size_t j = 0; j < this->actualBaseLearners; j++)
					value02 = value02 + theHypoWeight.get(j)*tmp.get(j);
				algorithmValue.set(i, value02);
			}
		}
	}
	
	void AlgorithmAdaBoost::classif(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmClass, sg::base::DataVector& algorithmValue)
	{
		eval(testData, algorithmValue);
		for (size_t i = 0; i < testData.getNrows(); i++)
		{
			algorithmClass.set(i, hValue(algorithmValue.get(i)));
		}
	}
	
	void AlgorithmAdaBoost::getAccuracyAndError(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, double* accruracy, double* error)
	{
		/* get the accuracy */
		size_t right = 0;
		sg::base::DataVector classOfAlgorithm(testDataClass.getSize());
		sg::base::DataVector valueOfAlgorithm(testDataClass.getSize());
		classOfAlgorithm.setAll(0.0);
		valueOfAlgorithm.setAll(0.0);
		classif(testData, classOfAlgorithm, valueOfAlgorithm);
		for(size_t i = 0; i < testData.getNrows(); i++)
		{
			if(classOfAlgorithm.get(i) == testDataClass.get(i))
				right = right + 1;
		}
		
		*accruracy = double(right)/double(classOfAlgorithm.getSize());
		
		/* get the interpolation error */
		valueOfAlgorithm.sub(testDataClass);
		valueOfAlgorithm.sqr();
		*error = valueOfAlgorithm.sum()/double(valueOfAlgorithm.getSize());
	}
	
	double AlgorithmAdaBoost::getAccuracy(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass)
	{
		size_t right = 0;
		sg::base::DataVector classOfAlgorithm(testDataClass.getSize());
		sg::base::DataVector valueOfAlgorithm(testDataClass.getSize());
		classif(testData, classOfAlgorithm, valueOfAlgorithm);
		for(size_t i = 0; i < testData.getNrows(); i++)
		{
			if(classOfAlgorithm.get(i) == testDataClass.get(i))
				right = right + 1;
		}
		return double(right)/double(classOfAlgorithm.getSize());
	}
	
	double AlgorithmAdaBoost::getError(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass)
	{
		sg::base::DataVector helper(testData.getNrows());
		eval(testData, helper);
		helper.sub(testDataClass);
		helper.sqr();
		return helper.sum()/double(helper.getSize());
	}
	
	void AlgorithmAdaBoost::evalBL(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmValue, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t yourBaseLearner)
	{
		if (yourBaseLearner > this->actualBaseLearners)
		{
			std::cout << "Please choose your base learner number not larger than " << this->actualBaseLearners << std::endl;
			abort();
		}
		else
		{
			sg::base::DataMatrix baseLearnerMatr(testData.getNrows(), yourBaseLearner);
			sg::base::DataVector p(this->dim);
			sg::base::OperationEval* operateEval = sg::op_factory::createOperationEval(*this->grid);
			sg::base::DataVector alphaT(storageAlpha.getNrows());
			sg::base::DataVector baseLearnerVec(testData.getNrows());
			for (size_t t = 0; t < yourBaseLearner; t++)
			{
				storageAlpha.getColumn(t, alphaT);
				baseLearnerVec.setAll(0.0);
				for (size_t i = 0; i < testData.getNrows(); i++)
				{
					testData.getRow(i, p);
					double value01 = operateEval->eval(alphaT, p);
					// when there is only one baselearner actually, we do as following, just use normal classify to get the value
					if (yourBaseLearner == 1)
						algorithmValue.set(i, value01);
					else
						baseLearnerVec.set(i, hValue(value01));
				}
				if (yourBaseLearner > 1)
					baseLearnerMatr.setColumn(t, baseLearnerVec);
			}
			if (yourBaseLearner > 1)
			{
				double value02;
				sg::base::DataVector tmp(yourBaseLearner);
				for (size_t i = 0; i < testData.getNrows(); i++)
				{
					value02 = 0.0;
					for(size_t k = 0; k < yourBaseLearner; k++)
						tmp.set(k, baseLearnerMatr.get(i, k));
					for(size_t j = 0; j < yourBaseLearner; j++)
						value02 = value02 + hypoWeight.get(j)*tmp.get(j);
					algorithmValue.set(i, value02);
				}
			}
		}
	}
	
	void AlgorithmAdaBoost::classifBL(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmClass, sg::base::DataVector& algorithmValue, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t yourBaseLearner)
	{
		evalBL(testData, algorithmValue, storageAlpha, hypoWeight, yourBaseLearner);
		for (size_t i = 0; i < testData.getNrows(); i++)
		{
			algorithmClass.set(i, hValue(algorithmValue.get(i)));
		}
	}
	
	void AlgorithmAdaBoost::getAccuracyAndErrorBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, double* accruracy, double* error, size_t yourBaseLearner)
	{
		/* get the accuracy */
		size_t right = 0;
		sg::base::DataVector classOfAlgorithm(testDataClass.getSize());
		sg::base::DataVector valueOfAlgorithm(testDataClass.getSize());
		classOfAlgorithm.setAll(0.0);
		valueOfAlgorithm.setAll(0.0);
		classifBL(testData, classOfAlgorithm, valueOfAlgorithm, storageAlpha, hypoWeight, yourBaseLearner);
		for(size_t i = 0; i < testData.getNrows(); i++)
		{
			if(classOfAlgorithm.get(i) == testDataClass.get(i))
				right = right + 1;
		}
		
		*accruracy = double(right)/double(testDataClass.getSize());
		
		/* get the interpolation error */
		valueOfAlgorithm.sub(testDataClass);
		valueOfAlgorithm.sqr();
		*error = valueOfAlgorithm.sum()/double(valueOfAlgorithm.getSize());
	}
	
	double AlgorithmAdaBoost::getAccuracyBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t yourBaseLearner)
	{
		size_t right = 0;
		sg::base::DataVector classOfAlgorithm(testDataClass.getSize());
		sg::base::DataVector valueOfAlgorithm(testDataClass.getSize());
		classOfAlgorithm.setAll(0.0);
		valueOfAlgorithm.setAll(0.0);
		classifBL(testData, classOfAlgorithm, valueOfAlgorithm, storageAlpha, hypoWeight, yourBaseLearner);
		for(size_t i = 0; i < testData.getNrows(); i++)
		{
			if(classOfAlgorithm.get(i) == testDataClass.get(i))
				right = right + 1;
		}
		return double(right)/double(testDataClass.getSize());
	}
	
	void AlgorithmAdaBoost::doRefinement(sg::base::DataVector& alpha_ada, sg::base::DataVector& rhs_ada, sg::base::DataVector& weight_ada, size_t ada_time)
	{
		sg::base::GridGenerator* myGenerator = this->grid->createGridGenerator();
	
		size_t refineNumber = 0.1 * (this->grid->getSize());
		sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&alpha_ada, refineNumber, 0.0);
		myGenerator->refine(myRefineFunc);
		delete myRefineFunc;

		sg::base::GridStorage* gridStorage_ada = this->grid->getStorage();
		this->gridPoint = gridStorage_ada->size();
		std::cout << std::endl;
		std::cout << "Refinement time step: " << ada_time << ", new grid size: " << this->gridPoint << ", refined number of grid points: " << refineNumber << std::endl;
		
		// extend alpha vector (new entries uninitialized)
		alpha_ada.resizeZero(this->gridPoint);
		rhs_ada.resizeZero(this->gridPoint);
		
		// calculate new alpha
		sg::datadriven::DMWeightMatrix WMatrix_ada(*grid, *data, *C, this->lamb, weight_ada);
		WMatrix_ada.generateb(*classes, rhs_ada);
		sg::solver::ConjugateGradients myCG_ada(this->imax, this->epsilon);
		myCG_ada.solve(WMatrix_ada, alpha_ada, rhs_ada, false, true, -1.0);
	}

	size_t AlgorithmAdaBoost::getActualBL()
	{
		return this->actualBaseLearners;
	}

}
}
