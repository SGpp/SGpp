#include <limits>

#include "LearnerOnlineSGD.hpp"

#include "base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp"
namespace sg
{

namespace datadriven
{

LearnerOnlineSGD::LearnerOnlineSGD(
    sg::datadriven::LearnerRegularizationType& regularization,
    const bool isRegression, const bool isVerbose) :
        Learner(regularization, isRegression, isVerbose),
        mainTrainDataset(NULL), mainClasses(NULL),
        testTrainDataset(NULL), testClasses(NULL),
        SGDCurrentIndex(0),
        numMainData(0), numMainDim(0),
        mainError(NULL),
        minibatchTrainDataset(NULL),
        minibatchClasses(NULL),
        minibatchError(NULL),
        errorPerBasisFunction(NULL),
        alphaAvg(NULL),
        currentGamma(0),
        currentRunIterations(0)

{}

void LearnerOnlineSGD::train(sg::base::DataMatrix& mainTrainDataset_,
                             sg::base::DataVector& mainClasses_,
                             sg::base::DataMatrix& testTrainDataset_,
                             sg::base::DataVector& testClasses_,
                             sg::base::RegularGridConfiguration& gridConfig,
                             sg::datadriven::LearnerOnlineSGDConfiguration& config_,
                             sg::base::AbstractRefinement& refinement
                            )
{
    /*
     * Initialization
     */

    using namespace sg::base;

    if (alpha_ != NULL)
        delete alpha_;

    if (grid_ != NULL)
        delete grid_;

    if (mainError != NULL)
        delete mainError;

    if (minibatchTrainDataset != NULL)
        delete minibatchTrainDataset;

    if (minibatchClasses != NULL)
        delete minibatchClasses;

    if (minibatchError != NULL)
        delete minibatchError;

    if (errorPerBasisFunction != NULL)
        delete errorPerBasisFunction;

    if (alphaAvg != NULL)
        delete alphaAvg;

    if (isTrained_ == true)
        isTrained_ = false;

    InitializeGrid(gridConfig);
    if (grid_ == NULL)
    {
        return;
    }

    /*
     * Members
     */

    mainTrainDataset = &mainTrainDataset_;
    mainClasses = &mainClasses_;

    testTrainDataset = &testTrainDataset_;
    testClasses = &testClasses_;

    config = config_;

    numMainData = mainTrainDataset->getNrows();
    numMainDim = mainTrainDataset->getNcols();

    mainError = new DataVector(numMainData);
    mainError->setAll(0.0);

    minibatchTrainDataset = new sg::base::DataMatrix(0, numMainDim);
    minibatchTrainDataset->addSize(config.minibatchSize);
    minibatchTrainDataset->setAll(0.0);

    minibatchClasses = new DataVector(config.minibatchSize);
    minibatchClasses->setAll(0.0);

    minibatchError = new DataVector(config.minibatchSize);
    minibatchError->setAll(0.0);

    errorPerBasisFunction = new DataVector(grid_->getSize());
    errorPerBasisFunction->setAll(0.0);

    alphaAvg = new DataVector(grid_->getSize());
    alphaAvg->setAll(0.0);

    currentGamma = config.gamma;

    /*
     * File descriptors
     */

    std::fstream ferr0, ferr1, ferr2, fgrid, fcoor;
    std::string dir = config.experimentDir;
    ferr0.open((dir + std::string("/ferr0")).c_str(),
               std::ios::out | std::ios::trunc);
    ferr1.open((dir + std::string("/ferr1")).c_str(),
               std::ios::out | std::ios::trunc);
    ferr2.open((dir + std::string("/ferr2")).c_str(),
               std::ios::out | std::ios::trunc);
    fgrid.open((dir + std::string("/fgrid")).c_str(),
               std::ios::out | std::ios::trunc);
    fcoor.open((dir + std::string("/fcoor")).c_str(),
               std::ios::out | std::ios::trunc);

    /*
     * SGD Order
     */

    for (size_t i = 0; i < numMainData; i++)
    {
        SGDIndexOrder.push_back(i);
    }

    std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
    SGDCurrentIndex = 0;

    /*
     * Refinement Functor
     */

    RefinementFunctor *functor = NULL;

    if(config.refinementType == "SURPLUS")
    {
        functor = new SurplusRefinementFunctor(alpha_,
                                               config.refinementNumPoints, 0.0);
    }

    if(config.refinementType == "MSE")
    {
        /*functor = new SurplusRefinementFunctor(alpha_,
        		        RefineConfig.refinementNumPoints, 0.0);*/
        /*functor = new SurplusRefinementFunctor(alpha_,
                RefineConfig.refinementNumPoints, 0.0);*/
        functor = new SurplusRefinementFunctor(errorPerBasisFunction,
                                               config.refinementNumPoints, 0.0);
    }

    if(config.refinementType == "WEIGHTED_ERROR_MINIBATCH")
    {
        functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
                  config.refinementNumPoints,
                  -std::numeric_limits<double>::infinity());
        WeightedErrorRefinementFunctor* wfunctor =
            (WeightedErrorRefinementFunctor*) functor;
        wfunctor->setTrainDataset(minibatchTrainDataset);
        wfunctor->setClasses(minibatchClasses);
        wfunctor->setErrors(minibatchError);
    }

    if(config.refinementType == "WEIGHTED_ERROR_ALL")
    {
        // FIXME: this case is not accounted for
        /*functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
         RefineConfig.refinementNumPoints, 0.0);
         WeightedErrorRefinementFunctor* wfunctor =
         (WeightedErrorRefinementFunctor*) functor;

         wfunctor->setTrainDataset(mainTrainDataset);
         wfunctor->setClasses(mainClasses);*/
    }

    if(config.refinementType == "PERSISTENT_ERROR")
    {
        functor = new PersistentErrorRefinementFunctor(alphaAvg, grid_,
                  config.refinementNumPoints,
                  -std::numeric_limits<double>::infinity());
        PersistentErrorRefinementFunctor* pfunctor =
            (PersistentErrorRefinementFunctor*) functor;
        pfunctor->setTrainDataset(minibatchTrainDataset);
        pfunctor->setClasses(minibatchClasses);
        pfunctor->setErrors(minibatchError);
    }

    if(config.refinementType == "CLASSIFICATION")
    {

        functor = new ClassificationRefinementFunctor(alpha_, grid_,
                  config.refinementNumPoints, 0.0);
        ClassificationRefinementFunctor* cfunctor =
            (ClassificationRefinementFunctor*) functor;
        cfunctor->setTrainDataset(mainTrainDataset);
        cfunctor->setClasses(mainClasses);
    }

    if(config.refinementType == "PREDICTIVE_ERROR_DIMENSION_REFINEMENT")
    {
       functor = new PredictiveRefinementDimensionIndicator(grid_, mainTrainDataset, mainError, config.refinementNumPoints, 0.0);
    }

    if (functor == NULL)
    {
        throw base::application_exception("Invalid refinement type");
    }


    /*
     * Perform runs
     */

    for (int countRun = 0; countRun < config.numRuns; countRun++)
    {

        currentRunIterations = 0;

        std::cout << "Run: " << countRun + 1 << std::endl;

        while (currentRunIterations < numMainData)
        {
            /*
             * Perform SGD
             */

            if(config.refinementCondition == "FIXED_NUMBER")
            {
                for (size_t i = 0;
                        i < config.numIterations;
                        i++)
                {
                    currentGamma = config.gamma*pow(1+config.gamma*config.lambda*(double)currentRunIterations, -2.0/3);
                    performSGDStep();
                    currentRunIterations++;
                }

            }

            if(config.refinementCondition == "SMOOTHED_ERROR_DECLINE")
            {

                double oldErrorSum = 0;
                double oldErrorLast = 0;

                double currentMinibatchError = 0;
                double ratio = 1;

                do
                {
                    // Run SGD
                    performSGDStep();
                    currentRunIterations++;

                    // Calculate smoothed error
                    currentMinibatchError = getError(minibatchTrainDataset, minibatchClasses, config.errorType, NULL);
                    if (smoothedErrorDeclineBuffer.size() >= config.smoothedErrorDeclineBufferSize)
                    {

                        // Calculate average of old minibatch errors
                        for (
                            std::list<double>::iterator it = smoothedErrorDeclineBuffer.begin();
                            it != smoothedErrorDeclineBuffer.end();
                            ++it
                        )
                        {
                            oldErrorSum += *it;
                        }

                        oldErrorLast = smoothedErrorDeclineBuffer.back();

                        // Update errorOnMinibatch
                        smoothedErrorDeclineBuffer.pop_back();
                        smoothedErrorDeclineBuffer.push_front(currentMinibatchError);

                        // Update ratio
                        ratio = (oldErrorLast - currentMinibatchError) / oldErrorSum;
                        ratio *= 1.0 / (double) config.smoothedErrorDeclineBufferSize;

                    }
                    else
                    {
                        smoothedErrorDeclineBuffer.push_front(currentMinibatchError);
                    }
                }
                while (ratio > config.smoothedErrorDecline);
            }

            /*
             * Error vectors
             */

            if (config.refinementType == "PERSISTENT_ERROR" ||
                    config.refinementType == "WEIGHTED_ERROR_MINIBATCH")
            {
                if (config.errorType == "ACCURACY")
                {
                    double accuracy = getError(minibatchTrainDataset,
                                               minibatchClasses,
                                               config.errorType,
                                               minibatchError);
                    if (accuracy == 1.0)
                    {
                        continue;
                    }
                }
            }

            if (config.refinementType == "PREDICTIVE_ERROR_DIMENSION_REFINEMENT")
            {
				(dynamic_cast<OnlinePredictiveRefinementDimension*>(&refinement))->free_refine(
						  grid_->getStorage(), static_cast<PredictiveRefinementDimensionIndicator*>(functor));
            }

            if (config.refinementType == "MSE")
            {
                getError(mainTrainDataset, mainClasses, config.errorType, mainError);
                mainError->sqr();
                OperationMultipleEval* eval = sg::op_factory::createOperationMultipleEval(*grid_,
                                              mainTrainDataset);
                eval->multTranspose(*mainError, *errorPerBasisFunction);
                delete eval;
            }

            /*
             * Output
             */

            double err0 = getError(minibatchTrainDataset, minibatchClasses, config.errorType, NULL);
            double err1 = getError(mainTrainDataset, mainClasses, config.errorType, NULL);
            double err2 = getError(testTrainDataset, testClasses, config.errorType, NULL);

            ferr0 << currentRunIterations << "," << err0 << std::endl;
            ferr1 << currentRunIterations << "," << err1 << std::endl;
            ferr2 << currentRunIterations << "," << err2 << std::endl;

            std::string grid_str;
            grid_->serialize(grid_str);
            fgrid << grid_str << std::endl;

            double percent = (double) currentRunIterations + countRun * (double) numMainData;
            percent /= (double) numMainData * config.numRuns;
            percent *= 100;
            if (percent > 100)
            {
                percent = 100;
            }

            printf("Accuracy: %2.2f%% (at %2.2f%%)\n", err1 * 100, percent);
            fflush(stdout);

            /*
             * Refinement
             */

            refinement.free_refine(grid_->getStorage(), functor);
            alpha_->resizeZero(grid_->getSize());
            alphaAvg->resizeZero(grid_->getSize());
        }
    }

    /*
     * Perform CG
     */

    std::cout << "Error before CG (" << config.errorType << "): " <<
    getError(mainTrainDataset, mainClasses, "ACCURACY", NULL) << std::endl;
    std::cout << "Error before CG (" << config.errorType << "): " <<
    getError(mainTrainDataset, mainClasses, "MSE", NULL) << std::endl;

    sg::solver::ConjugateGradients *cg = new sg::solver::ConjugateGradients(
                                             config.CG_max, config.CG_eps);

    sg::base::OperationMatrix *C_ = sg::op_factory::createOperationIdentity(
                                        *this->grid_);
    sg::datadriven::DMSystemMatrix matrix(*grid_, *mainTrainDataset, *C_,
                                          config.lambda);

    sg::base::DataVector b(alpha_->getSize());
    matrix.generateb(*mainClasses, b);

    cg->solve(matrix, *alpha_, b, true, false);

    std::cout << "Error after CG (" << config.errorType << "): " <<
    getError(mainTrainDataset, mainClasses, "ACCURACY", NULL) << std::endl;
    std::cout << "Error after CG (" << config.errorType << "): " <<
    getError(mainTrainDataset, mainClasses, "MSE", NULL) << std::endl;

    std::cout << "Error on test (" << config.errorType << "): " <<
    getError(testTrainDataset, mainClasses, "ACCURACY", NULL) << std::endl;
    std::cout << "Error on test (" << config.errorType << "): " <<
    getError(testTrainDataset, mainClasses, "MSE", NULL) << std::endl;

    /*
     * Clean up
     */

    isTrained_ = true;

    ferr0.close();
    ferr1.close();
    ferr2.close();
    fgrid.close();
    fcoor.close();

    delete C_;
    delete cg;
    delete functor;
}

double LearnerOnlineSGD::getError(sg::base::DataMatrix* trainDataset,
                                  sg::base::DataVector* classes, std::string errorType, sg::base::DataVector* error = NULL)
{
    using namespace sg::base;

    size_t numData = trainDataset->getNrows();

    bool cleanup = false;
    if( error == NULL )
    {
        error = new DataVector(numData);
        cleanup = true;
    }

    DataVector result(numData);
    OperationMultipleEval* eval = sg::op_factory::createOperationMultipleEval(*grid_, trainDataset);

    eval->mult(*alphaAvg, result);
    //eval->mult(*alpha_, result);
    delete eval;

    double res = -1.0;

    if(errorType == "MSE")
    {
        for (unsigned int i = 0; i < numData; i++)
        {
            error->set(i, classes->get(i) - result.get(i));
        }

        // Error
        double sum = 0;
        for (unsigned int i = 0; i < numData; i++)
        {
            sum += error->get(i) * error->get(i);
        }

        res = (sum / (double) numData);

    }
    if(errorType == "ACCURACY")
    {
        unsigned int correct = 0;
        for (unsigned int i = 0; i < numData; i++)
        {
            correct += (result.get(i) < 0) == (classes->get(i) < 0) ? 1 : 0;
            error->set(i, classes->get(i) - result.get(i));

        }
        res = static_cast<double>(correct) / static_cast<double>(numData);

    }

    if (cleanup)
    {
        delete error;
    }

    return res;
}

void LearnerOnlineSGD::performSGDStep()
{
    using namespace sg::base;

    size_t numCoeff = grid_->getStorage()->size();

    // Get x and y pair
    DataVector x(numMainDim);
    mainTrainDataset->getRow(SGDIndexOrder[SGDCurrentIndex], x);
    double y = mainClasses->get(SGDIndexOrder[SGDCurrentIndex]);

    // Store in minibatch
    pushMinibatch(x, y);


	// Update SGDCurrentIndex
	if (SGDCurrentIndex == SGDIndexOrder.size() - 1) {
		std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
		SGDCurrentIndex = 0;
	} else {
		SGDCurrentIndex++;
	}

    // Calculate delta^n according to [Maier BA, 5.10]:
    // b_k^T * alpha^n - y_k
    double tmp1 = grid_->eval(*alpha_, x) - y;

    // delta^n = 2 * gamma * (b_k * tmp1 + lambda * a^n)
    DataVector delta(*alpha_);
    DataVector unit_alpha(numCoeff);
    unit_alpha.setAll(0.0);

    DataVector singleAlpha(1);
    singleAlpha[0] = 1.0;

    DataMatrix dm(x.getPointer(), 1, x.getSize());
    OperationMultipleEval *multEval = sg::op_factory::createOperationMultipleEval(*grid_, &dm);
    multEval->multTranspose(singleAlpha, delta);
    delete multEval;
    alpha_->mult(1-currentGamma * config.lambda);
    alpha_->axpy(-currentGamma * tmp1, delta);

    //double mu = 0.1; // boring exponential smoothing

    // L. Bottou exciting smoothing
    size_t t1 = (currentRunIterations > numMainDim+1) ? currentRunIterations - numMainDim : 1;
    size_t t2 = (currentRunIterations > numMainData+1) ? currentRunIterations - numMainData: 1;
    double mu = (t1>t2) ? static_cast<double>(t1) : static_cast<double>(t2);
    mu = 1.0/mu;

    alphaAvg->mult(1-mu);
    alphaAvg->axpy(mu, *alpha_);

    /*for (size_t i = 0; i < numCoeff; i++) {
    	unit_alpha[i] = 1;
    	delta[i] = grid_->eval(unit_alpha, x) * tmp1;
    	delta[i] += lambda * (*alpha_)[i];
    	delta[i] *= 2 * gamma;
    	unit_alpha[i] = 0;
}

    // update alpha
    // a^{n+1} = a^n - delta^n
    for (size_t i = 0; i < numCoeff; i++) {
    	(*alpha_)[i] = (*alpha_)[i] - delta[i];
}*/
}

void LearnerOnlineSGD::pushMinibatch(sg::base::DataVector& x, double y)
{
    static size_t next_idx = 0;
    if (minibatchTrainDataset->getUnused() > 0)
    {
        minibatchTrainDataset->appendRow(x);
        (*minibatchClasses)[next_idx] = y;
    }
    else
    {
        minibatchTrainDataset->setRow(next_idx, x);
        (*minibatchClasses)[next_idx] = y;
    }
    next_idx = (next_idx + 1) % config.minibatchSize;
}


LearnerOnlineSGD::~LearnerOnlineSGD()
{
    if (mainError != NULL)
        delete mainError;

    if (minibatchTrainDataset != NULL)
        delete minibatchTrainDataset;

    if (minibatchClasses != NULL)
        delete minibatchClasses;

    if (minibatchError != NULL)
        delete minibatchError;

    if (errorPerBasisFunction != NULL)
        delete errorPerBasisFunction;

    if (alphaAvg != NULL)
        delete alphaAvg;
}

}

}
