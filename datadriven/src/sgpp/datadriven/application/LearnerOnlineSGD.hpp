/* ****************************************************************************
 * Copyright (C) 2014 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Maxim Schmidt (maxim.schmidt@tum.de)
#ifndef LEARNERONLINESGD_HPP
#define LEARNERONLINESGD_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

//TODO remove meta header usage
#include "sgpp_base.hpp"
#include "sgpp_datadriven.hpp"
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimension.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimensionOld.hpp>
#include <sgpp/parallel/tools/TypesParallel.hpp>


namespace sg
{

namespace datadriven
{


struct LearnerOnlineSGDConfiguration
{
    std::string refinementCondition;
    size_t numIterations;
    size_t smoothedErrorDeclineBufferSize;

    std::string refinementType;
    int refinementNumPoints;

    double lambda;
    double gamma;

    std::string errorType;

    int numRuns;

    std::string experimentDir;

    size_t minibatchSize;

    int CG_max;
    double CG_eps;
    double smoothedErrorDecline;

    std::string hashRefinementType;
};

class LearnerOnlineSGD: public sg::datadriven::Learner
{

public:
    LearnerOnlineSGD(sg::datadriven::LearnerRegularizationType& regularization,
                     const bool isRegression, const bool isVerbose = true);



    virtual void train(sg::base::DataMatrix& mainTrainDataset_,
                       sg::base::DataVector& mainClasses_,

                       sg::base::DataMatrix& testTrainDataset_,
                       sg::base::DataVector& testClasses_,

                       sg::base::RegularGridConfiguration& gridConfig,
                       sg::datadriven::LearnerOnlineSGDConfiguration& config
                       );

    virtual ~LearnerOnlineSGD();
protected:
    static const sg::parallel::VectorizationType vecType_;

private:
    sg::base::DataMatrix* mainTrainDataset;
    sg::base::DataVector* mainClasses;

    sg::base::DataMatrix* testTrainDataset;
    sg::base::DataVector* testClasses;

    sg::datadriven::LearnerOnlineSGDConfiguration config;

    size_t SGDCurrentIndex;
    std::vector<size_t> SGDIndexOrder;

    size_t numMainData;
    size_t numMainDim;

    sg::base::HashRefinement* hash_refinement;
    sg::base::OnlinePredictiveRefinementDimension* online_refinement;
    sg::base::OnlinePredictiveRefinementDimensionOld* online_refinement_old;


    double getError(sg::base::DataMatrix* trainDataset,
                    sg::base::DataVector* classes,
                    std::string errorType,
                    sg::base::DataVector* error, bool useEvalVectorized);

    void pushMinibatch(sg::base::DataVector& x, double y);
    void performSGDStep();

    /*
     * Other
     */

    sg::base::DataVector* mainError;

    sg::base::DataMatrix* minibatchTrainDataset;
    sg::base::DataVector* minibatchClasses;
    sg::base::DataVector* minibatchError;

    sg::base::DataVector* errorPerBasisFunction;
    sg::base::DataVector* alphaAvg;

    std::list<double> smoothedErrorDeclineBuffer;
    double currentGamma;
    size_t currentRunIterations;
};

}
}

#endif /* LEARNERSGD_HPP */
