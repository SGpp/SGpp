// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ALGORITHMADABOOSTVECTORIZEDIDENTITY_HPP
#define ALGORITHMADABOOSTVECTORIZEDIDENTITY_HPP

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp>

#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    /*
     * Algorithm for Adaboosting
     * This algorithm is to train Train base learner according to sample distribution and obtain hypothesis
     * get the hypothesis weight
     * Then combine hypothesis linearly
     *
     * The main idea behind the algorithm is to get a better accuracy in classify dateset according to the training dataset
     *
     * Vectorized MultipleEvalulation Operations are used!
     */
    class AlgorithmAdaBoostVectorizedIdentity : public SGPP::datadriven::AlgorithmAdaBoostBase {
      protected:
        /// Vectorization mode, possible values are SSE, AVX, OCL, ArBB
        VectorizationType vecMode;

        virtual void alphaSolver(double& lambda, SGPP::base::DataVector& weight, SGPP::base::DataVector& alpha, bool final);

      public:
        /**
         * Std-Constructor
         *
         * @param SparseGrid reference to the sparse grid
         * @param gridType reference to the of grid type(1 = Linear Grid, 2 = LinearL0Boundary Grid, 3 = ModLinear Grid)
         * @param gridLevel reference to the level of grid
         * @param trainData reference to the training dataset
         * @param trainDataClass reference to the class of training dataset
         * @param NUM the number of baselearner for Adaboosting
         * @param lambda the regularisation parameter
         * @param IMAX the parameter for ConjugateGradients
         * @param eps the parameter for ConjugateGradients
         * @param IMAX_final the parameter for ConjugateGradients used for last refinement step
         * @param eps_final the parameter for ConjugateGradients used for last refinement step
         * @param firstLabel one label from training dataset
         * @param secondLabel another label from training dataset
         * @param threshold the parameter for predicting a class
         * @param maxLambda the max lambda used in searching optimal lambda
         * @param minLambda the min lambda used in searching optimal lambda
         * @param searchNum the searching times used in searching for optimal lambda
         * @param refine the judgement of refine
         * @param refineMode Select the refine mode
         * @param refineNum the Number of refinement with a certain percentage of Grid points
         * @param numberOfAda the number of Grid points to refine
         * @param percentOfAda the percentage of Grid points to refine
         * @param vecMode vectorization mode, possible values are SSE, AVX, OCL, ArBB
         * @param mode the adaboost type to choose
         */
        AlgorithmAdaBoostVectorizedIdentity(SGPP::base::Grid& SparseGrid, size_t gridType, int gridLevel, SGPP::base::DataMatrix& trainData, SGPP::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, size_t IMAX_final, double eps_final, double firstLabel, double secondLabel, double threshold, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, double percentOfAda, VectorizationType vecMode, size_t mode);

        /**
         * Std-Deconstructor
         */
        virtual ~AlgorithmAdaBoostVectorizedIdentity();
    };

  }

}

#endif /* ALGORITHMADABOOSTVECTORIZEDIDENTITY_HPP */