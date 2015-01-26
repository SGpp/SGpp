/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (songz@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)
// @auther Alexander Heinecke (alexander.heinecke@mytum.de)

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostIdentity.hpp>
#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    AlgorithmAdaBoostIdentity::AlgorithmAdaBoostIdentity(SGPP::base::Grid& SparseGrid, size_t gridType, SGPP::base::HashGenerator::level_t gridLevel, SGPP::base::DataMatrix& trainData, SGPP::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, size_t IMAX_final, double eps_final, double firstLabel, double secondLabel, double threshold, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, double percentOfAda, size_t mode) : AlgorithmAdaBoostBase(SparseGrid, gridType, gridLevel, trainData, trainDataClass, NUM, lambda, IMAX, eps, IMAX_final, eps_final, firstLabel, secondLabel, threshold, maxLambda, minLambda, searchNum, refine, refineMode, refineNum, numberOfAda, percentOfAda, mode) {
    }

    AlgorithmAdaBoostIdentity::~AlgorithmAdaBoostIdentity() {
    }

    void AlgorithmAdaBoostIdentity::alphaSolver(double& lambda, SGPP::base::DataVector& weight, SGPP::base::DataVector& alpha, bool final) {
      SGPP::base::OperationMatrix* C = SGPP::op_factory::createOperationIdentity(*this->grid);
      SGPP::datadriven::DMWeightMatrix WMatrix(*this->grid, *this->data, *C, lambda, weight);
      SGPP::base::DataVector rhs(alpha.getSize());
      WMatrix.generateb(*this->classes, rhs);

      if (final) {
        SGPP::solver::ConjugateGradients myCG(this->imax_final, this->epsilon_final);
        myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
      } else {
        SGPP::solver::ConjugateGradients myCG(this->imax, this->epsilon);
        myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
      }

      delete C;
    }

  }
}
