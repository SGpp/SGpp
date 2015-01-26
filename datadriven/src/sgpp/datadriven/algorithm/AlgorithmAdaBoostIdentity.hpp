/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Zhongwen Song (songz@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)
// @auther Alexander Heinecke (alexander.heinecke@mytum.de)

#ifndef ALGORITHMADABOOSTIDENTITY_HPP
#define ALGORITHMADABOOSTIDENTITY_HPP

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    /*
     * Algorithm for Adaboosting
     *
     * For training the grid, a sparse grid standard approach is taken.
     */
    class AlgorithmAdaBoostIdentity : public AlgorithmAdaBoostBase {
      protected:
        /**
         * Performs a solver to get alpha use DMWeightMatrix as the System Matrix
        *
         * @param lambda the regularisation parameter
        * @param weight the weights of examples
        * @param alpha output the coefficients of the sparse grid's basis functions
        * @param final judgement the final step of this base learner
        */
        virtual void alphaSolver(double& lambda, SGPP::base::DataVector& weight, SGPP::base::DataVector& alpha, bool final);

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
        * @param mode the adaboost type to choose
         */
        AlgorithmAdaBoostIdentity(SGPP::base::Grid& SparseGrid, size_t gridType, SGPP::base::HashGenerator::level_t gridLevel, SGPP::base::DataMatrix& trainData, SGPP::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, size_t IMAX_final, double eps_final, double firstLabel, double secondLabel, double threshold, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, double percentOfAda, size_t mode);

        /**
         * Std-Deconstructor
         */
        virtual ~AlgorithmAdaBoostIdentity();
    };

  }
}
#endif



