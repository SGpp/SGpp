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
    }
    
    AlgorithmAdaBoost::~AlgorithmAdaBoost()
    {
        delete this->C;
    }

	double AlgorithmAdaBoost::hValue(double realValue)
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
	}
    
    void AlgorithmAdaBoost::doAdaBoost(sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t* actualBN)
    {
        *actualBN = 0;
		sg::base::DataVector weight(this->numData);
        for (size_t i = 0; i < this->numData;i++)
        {
			weight.set(i, 1.0/this->numData);
        }
		sg::base::GridStorage* gridStorage = this->grid->getStorage();
            // create coefficient vector
		sg::base::DataVector alpha(gridStorage->size());
            // create weight error vector(ε)
		sg::base::DataVector weighterror(this->numBaseLearners);
		sg::base::DataVector rhs(gridStorage->size());
		sg::base::OperationEval* opEval = sg::GridOperationFactory::createOperationEval(*this->grid);
            // to store certain train data point
		sg::base::DataVector p(this->dim);
            // create vector to store the hypothesis of the training data according to certain alpha vector(base learner)
		sg::base::DataVector newclasses(this->numData);
            // create vector to store the identity of comparing hypothesis of training data with class of
		sg::base::DataVector identity(this->numData);

		sg::base::DataVector tmpweight(this->numData);
        for (size_t count = 0; count < this->numBaseLearners; count++){
			(*actualBN)++;
            alpha.setAll(0.0);
            
			sg::datadriven::DMWeightMatrix WMatrix(*grid, *data, *C, this->lamb, weight);
            WMatrix.generateb(*classes, rhs);

			sg::solver::ConjugateGradients myCG(this->imax, this->epsilon);
            myCG.solve(WMatrix, alpha, rhs, false, true, -1.0);
            
            storageAlpha.setColumn(count, alpha);

			for (size_t i = 0; i < this->numData; i++){
                this->data->getRow(i, p);
                double value = opEval->eval(alpha, p);
                newclasses.set(i, hValue(value));
            }
            
            for (size_t i = 0; i < this->numData; i++){
                if (newclasses.get(i) == this->classes->get(i)){
                    identity.set(i, 0.0);
                }else{
                    identity.set(i, 1.0);
                }
            }

			// to judge the classification
            if (identity.sum() == 0){
                std::cout << std::endl << "There is no wrong classification any more after " << count + 1 << " iterations" << std::endl;
                break;
            }

            weighterror.set(count, weight.dotProduct(identity));
            hypoWeight.set(count, 0.5*(log((1.0-weighterror.get(count))/weighterror.get(count))));

			std::cout << count << " " << hypoWeight.get(count)  << std::endl; 

            for (size_t i = 0; i < this->numData; i++){
                double tmp = weight.get(i)*exp(-hypoWeight.get(count) * newclasses.get(i)*this->classes->get(i));
                tmpweight.set(i, tmp);
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
		storageOfAlpha.setAll(0.0);
		theHypoWeight.setAll(0.0);
		size_t actualBN;
		doAdaBoost(storageOfAlpha, theHypoWeight, &actualBN);
		sg::base::DataMatrix baseLearnerMatr(testData.getNrows(), actualBN);
		sg::base::DataVector finalHypoClasses(actualBN);
		sg::base::DataVector p(this->dim);
		sg::base::OperationEval* operateEval = sg::GridOperationFactory::createOperationEval(*this->grid);
		sg::base::DataVector alphaT(storageOfAlpha.getNrows());
		sg::base::DataVector baseLearnerVec(testData.getNrows());
        for (size_t t = 0; t < actualBN; t++){
			storageOfAlpha.getColumn(t, alphaT);
			baseLearnerVec.setAll(0.0);
            for (size_t i = 0; i < testData.getNrows(); i++){
                testData.getRow(i, p);
                double value = operateEval->eval(alphaT, p);
				// when there is only one baselearner actually, we do as following, just use normal classify to get the value
				if (actualBN == 1)
					algorithmValue.set(i, value);
				else
					baseLearnerVec.set(i, hValue(value));
            }
			if (actualBN > 1)
				baseLearnerMatr.setColumn(t, baseLearnerVec);
        }
		if (actualBN > 1)
		{
			double value;
			sg::base::DataVector tmp(actualBN);
			for (size_t i = 0; i < testData.getNrows(); i++)
			{
				value = 0.0;
				for(size_t k = 0; k < actualBN; k++)
					tmp.set(k, baseLearnerMatr.get(i, k));
				for(size_t j = 0; j < actualBN; j++)
					value = value + theHypoWeight.get(j)*tmp.get(j);
				algorithmValue.set(i, value);
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
        double right = 0.0;
		sg::base::DataVector classOfAlgorithm(testDataClass.getSize());
		sg::base::DataVector valueOfAlgorithm(testDataClass.getSize());
		classOfAlgorithm.setAll(0.0);
		valueOfAlgorithm.setAll(0.0);
        classif(testData, classOfAlgorithm, valueOfAlgorithm);
        for(size_t i = 0; i < testData.getNrows(); i++)
        {
            if(classOfAlgorithm.get(i) == testDataClass.get(i))
                right = right + 1.0;
        }

        *accruracy = right/classOfAlgorithm.getSize();

		/* get the interpolation error */
        valueOfAlgorithm.sub(testDataClass);
        valueOfAlgorithm.sqr();
        *error = valueOfAlgorithm.sum()/valueOfAlgorithm.getSize();
    }

    double AlgorithmAdaBoost::getAccuracy(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass)
    {
        double right = 0.0;
		sg::base::DataVector classOfAlgorithm(testDataClass.getSize());
		sg::base::DataVector valueOfAlgorithm(testDataClass.getSize());
        classif(testData, classOfAlgorithm, valueOfAlgorithm);
        for(size_t i = 0; i < testData.getNrows(); i++)
        {
            if(classOfAlgorithm.get(i) == testDataClass.get(i))
                right = right + 1.0;
        }
        return right/classOfAlgorithm.getSize();
    }

    double AlgorithmAdaBoost::getError(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass)
    {
		sg::base::DataVector helper(testData.getNrows());
        eval(testData, helper);
        helper.sub(testDataClass);
        helper.sqr();
        return helper.sum()/helper.getSize();
    }

    void AlgorithmAdaBoost::evalBL(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmValue, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t actualBN, size_t yourBaseLearner)
    {
		if (yourBaseLearner > actualBN)
		{
			std::cout << "Please choose your base learner number not larger than " << actualBN << endl;
			abort();
		}
		else
		{
			sg::base::DataMatrix baseLearnerMatr(testData.getNrows(), yourBaseLearner);
			sg::base::DataVector finalHypoClasses(yourBaseLearner);
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
					double value = operateEval->eval(alphaT, p);
					// when there is only one baselearner actually, we do as following, just use normal classify to get the value
					if (yourBaseLearner == 1)
						algorithmValue.set(i, value);
					else
						baseLearnerVec.set(i, hValue(value));
				}
				if (yourBaseLearner > 1)
					baseLearnerMatr.setColumn(t, baseLearnerVec);
			}
			if (yourBaseLearner > 1)
			{
				double value;
				sg::base::DataVector tmp(yourBaseLearner);
				for (size_t i = 0; i < testData.getNrows(); i++)
				{
					value = 0.0;
					for(size_t k = 0; k < yourBaseLearner; k++)
						tmp.set(k, baseLearnerMatr.get(i, k));
					for(size_t j = 0; j < yourBaseLearner; j++)
						value = value + hypoWeight.get(j)*tmp.get(j);
					algorithmValue.set(i, value);
				}
			}
		}
	}

    void AlgorithmAdaBoost::classifBL(sg::base::DataMatrix& testData, sg::base::DataVector& algorithmClass, sg::base::DataVector& algorithmValue, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t actualBN, size_t yourBaseLearner)
    {
        evalBL(testData, algorithmValue, storageAlpha, hypoWeight, actualBN, yourBaseLearner);
        for (size_t i = 0; i < testData.getNrows(); i++)
        {
            algorithmClass.set(i, hValue(algorithmValue.get(i)));
        }
    }

    void AlgorithmAdaBoost::getAccuracyAndErrorBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, double* accruracy, double* error, size_t actualBN, size_t yourBaseLearner)
    {
		/* get the accuracy */
        double right = 0.0;
		sg::base::DataVector classOfAlgorithm(testDataClass.getSize());
		sg::base::DataVector valueOfAlgorithm(testDataClass.getSize());
		classOfAlgorithm.setAll(0.0);
		valueOfAlgorithm.setAll(0.0);
        classifBL(testData, classOfAlgorithm, valueOfAlgorithm, storageAlpha, hypoWeight, actualBN, yourBaseLearner);
        for(size_t i = 0; i < testData.getNrows(); i++)
        {
            if(classOfAlgorithm.get(i) == testDataClass.get(i))
                right = right + 1.0;
        }

        *accruracy = right/testDataClass.getSize();

		/* get the interpolation error */
        valueOfAlgorithm.sub(testDataClass);
        valueOfAlgorithm.sqr();
        *error = valueOfAlgorithm.sum()/valueOfAlgorithm.getSize();
    }

    double AlgorithmAdaBoost::getAccuracyBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t actualBN, size_t yourBaseLearner)
    {
        double right = 0.0;
		sg::base::DataVector classOfAlgorithm(testDataClass.getSize());
		sg::base::DataVector valueOfAlgorithm(testDataClass.getSize());
		classOfAlgorithm.setAll(0.0);
		valueOfAlgorithm.setAll(0.0);
        classifBL(testData, classOfAlgorithm, valueOfAlgorithm, storageAlpha, hypoWeight, actualBN, yourBaseLearner);
        for(size_t i = 0; i < testData.getNrows(); i++)
        {
            if(classOfAlgorithm.get(i) == testDataClass.get(i))
                right = right + 1.0;
        }
        return right/testDataClass.getSize();
    }

    double AlgorithmAdaBoost::getErrorBL(sg::base::DataMatrix& testData, sg::base::DataVector& testDataClass, sg::base::DataMatrix& storageAlpha, sg::base::DataVector& hypoWeight, size_t actualBN, size_t yourBaseLearner)
    {
		sg::base::DataVector helper(testData.getNrows());
        evalBL(testData, helper, storageAlpha, hypoWeight, actualBN, yourBaseLearner);
        helper.sub(testDataClass);
        helper.sqr();
        return helper.sum()/helper.getSize();
    }
}
}
