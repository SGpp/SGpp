/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (tjsongzw@gmail.com)
// @author Benjamin (pehersto@in.tum.de)
#include "algorithm/datadriven/AlgorithmAdaBoost.hpp"
#include "exception/algorithm_exception.hpp"
#include "basis/operations_factory.hpp"

namespace sg
{
namespace datadriven
{
	AlgorithmAdaBoost::AlgorithmAdaBoost(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, double firstLabel, double secondLabel)
    {
		sg::base::GridStorage* gridStorage = SparseGrid.getStorage();
        this->grid = &SparseGrid;
        this->lamb = lambda;
        this->data = &trainData;
        this->classes = &trainDataClass;
        this->numData = trainData.getNrows();
        this->dim = gridStorage->dim();
        this->numBaseLearners = NUM;
		this->C = sg::GridOperationFactory::createOperationIdentity(SparseGrid);
        this->imax = IMAX;
        this->epsilon = eps;
		this->labelOne = firstLabel;
		this->labelTwo = secondLabel;
		this->actualBaseLearners = 0;
    }
    
    AlgorithmAdaBoost::~AlgorithmAdaBoost()
    {
        delete this->C;
    }

    void AlgorithmAdaBoost::doAdaBoost(sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, sg::base::DataMatrix& weights, sg::base::DataMatrix& decision)
    {
		sg::base::DataVector weight(this->numData);
		weight.setAll(1.0/double(this->numData));

		sg::base::GridStorage* gridStorage = this->grid->getStorage();
            // create coefficient vector
		sg::base::DataVector alpha(gridStorage->size());
            // create weight error vector(ε)
		sg::base::DataVector weighterror(this->numBaseLearners);
		weighterror.setAll(0.0);
		sg::base::DataVector rhs(gridStorage->size());
		sg::base::OperationEval* opEval = sg::GridOperationFactory::createOperationEval(*this->grid);
            // to store certain train data point
		sg::base::DataVector p(this->dim);
            // create vector to store the hypothesis of the training data according to certain alpha vector(base learner)
		sg::base::DataVector newclasses(this->numData);
            // create vector to store the identity of comparing hypothesis of training data with class of
		sg::base::DataVector identity(this->numData);

		sg::base::DataVector tmpweight(this->numData);
		    // stop criterion
		// bool stop;
		sg::base::DataMatrix baseLearnerMatr(this->numData, this->numBaseLearners);
		baseLearnerMatr.setAll(0.0);

        for (size_t count = 0; count < this->numBaseLearners; count++){
			// stop = true;
			(this->actualBaseLearners)++;
			weights.setColumn(count, weight);
			alpha.setAll(0.0);
			sg::datadriven::DMWeightMatrix WMatrix(*grid, *data, *C, this->lamb, weight);
            WMatrix.generateb(*classes, rhs);

			sg::solver::ConjugateGradients myCG(this->imax, this->epsilon);
            myCG.solve(WMatrix, alpha, rhs, false, true, -1.0);
            
            storageAlpha.setColumn(count, alpha);

			for (size_t i = 0; i < this->numData; i++){
                this->data->getRow(i, p);
                double value01 = opEval->eval(alpha, p);
                newclasses.set(i, hValue(value01));
            }
			baseLearnerMatr.setColumn(count, newclasses);
            
            for (size_t i = 0; i < this->numData; i++){
                if (newclasses.get(i) == this->classes->get(i)){
                    identity.set(i, 0.0);
					decision.set(i, count, 1.0);
                }else{
                    identity.set(i, 1.0);
					decision.set(i, count, 0.0);
                }
            }
			// calculate the weighted train error rate
            weighterror.set(count, weight.dotProduct(identity));
			// to judge the classification
			if (weighterror.get(count) >= 0.5)
			{
                std::cout << std::endl << "The training error rate exceeds 0.5 after " << count + 1 << " iterations" << std::endl;
				(this->actualBaseLearners)--;
				break;
			}
			if (identity.sum() == 0)
			{
                std::cout << std::endl << "The training error rate is 0 after  " << count + 1 << " iterations" << std::endl;
				(this->actualBaseLearners)--;
				break;
			}

			// calculate the weight of this weak classify
            hypoWeight.set(count, 0.5*(log((1.0-weighterror.get(count))/weighterror.get(count))));

			/*
			if (this->actualBaseLearners == 1)
				for (size_t i = 0; i < this->numData; i++)
				{
					if (baseLearnerMatr.get(i, 0) != this->classes->get(i))
					{
						stop = false;
						break;
					}
				}
			else
			{
				double value02;
				sg::base::DataVector tmp(this->actualBaseLearners);
				for (size_t i = 0; i < this->numData; i++)
				{
					value02 = 0.0;
					for(size_t k = 0; k < this->actualBaseLearners; k++)
						tmp.set(k, baseLearnerMatr.get(i, k));
					for(size_t j = 0; j < this->actualBaseLearners; j++)
						value02 = value02 + hypoWeight.get(j)*tmp.get(j);
					if(hValue(value02) != this->classes->get(i))
					{
						stop = false;
						break;
					}
				}
			}
            if (stop)
			{
                std::cout << std::endl << "There is no wrong classification any more after " << count + 1 << " iterations" << std::endl;
                break;
			}
			*/

			double helper;
            for (size_t i = 0; i < this->numData; i++){
                helper = weight.get(i)*exp(-hypoWeight.get(count) * newclasses.get(i)*this->classes->get(i));
                tmpweight.set(i, helper);
            }

            double sum = tmpweight.sum();
            tmpweight.mult(1.0/sum);
                // get new weights vector
            weight = tmpweight;
        }
    }

    void AlgorithmAdaBoost::eval(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmValue)
    {
		sg::base::GridStorage* gridStorage = this->grid->getStorage();
		sg::base::DataMatrix storageOfAlpha(gridStorage->size(), this->numBaseLearners);
		sg::base::DataVector theHypoWeight(this->numBaseLearners);
		sg::base::DataMatrix weightsMatrix(this->numData, this->numBaseLearners);
		sg::base::DataMatrix decisionMatrix(this->numData, this->numBaseLearners);
		storageOfAlpha.setAll(0.0);
		theHypoWeight.setAll(0.0);
		weightsMatrix.setAll(0.0);
		decisionMatrix.setAll(0.0);
		doAdaBoost(storageOfAlpha, theHypoWeight, weightsMatrix, decisionMatrix);
		sg::base::DataMatrix baseLearnerMatr(testData.getNrows(), this->actualBaseLearners);
		sg::base::DataVector p(this->dim);
		sg::base::OperationEval* operateEval = sg::GridOperationFactory::createOperationEval(*this->grid);
		sg::base::DataVector alphaT(storageOfAlpha.getNrows());
		sg::base::DataVector baseLearnerVec(testData.getNrows());
        for (size_t t = 0; t < this->actualBaseLearners; t++){
			storageOfAlpha.getColumn(t, alphaT);
			baseLearnerVec.setAll(0.0);
            for (size_t i = 0; i < testData.getNrows(); i++){
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
			std::cout << "Please choose your base learner number not larger than " << this->actualBaseLearners << endl;
			abort();
		}
		else
		{
			sg::base::DataMatrix baseLearnerMatr(testData.getNrows(), yourBaseLearner);
			sg::base::DataVector p(this->dim);
			sg::base::OperationEval* operateEval = sg::GridOperationFactory::createOperationEval(*this->grid);
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

    double AlgorithmAdaBoost::getErrorBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t yourBaseLearner)
    {
		sg::base::DataVector helper(testData.getNrows());
        evalBL(testData, helper, storageAlpha, hypoWeight, yourBaseLearner);
        helper.sub(testDataClass);
        helper.sqr();
        return helper.sum()/double(helper.getSize());
    }

	size_t AlgorithmAdaBoost::getActualBL()
	{
		return this->actualBaseLearners;
	}
}
}
