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

namespace sg {
  namespace datadriven {

    AlgorithmAdaBoostIdentity::AlgorithmAdaBoostIdentity(sg::base::Grid& SparseGrid, size_t gridType, sg::base::HashGenerator::level_t gridLevel, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, size_t IMAX_final, double eps_final, double firstLabel, double secondLabel, double threshold, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, double percentOfAda, size_t mode) : AlgorithmAdaBoostBase(SparseGrid, gridType, gridLevel, trainData, trainDataClass, NUM, lambda, IMAX, eps, IMAX_final, eps_final, firstLabel, secondLabel, threshold, maxLambda, minLambda, searchNum, refine, refineMode, refineNum, numberOfAda, percentOfAda, mode) {
    }

    AlgorithmAdaBoostIdentity::~AlgorithmAdaBoostIdentity() {
    }

    void AlgorithmAdaBoostIdentity::alphaSolver(double& lambda, sg::base::DataVector& weight, sg::base::DataVector& alpha, bool final) {
      sg::base::OperationMatrix* C = sg::op_factory::createOperationIdentity(*this->grid);
      sg::datadriven::DMWeightMatrix WMatrix(*this->grid, *this->data, *C, lambda, weight);
      sg::base::DataVector rhs(alpha.getSize());
      WMatrix.generateb(*this->classes, rhs);

      if (final) {
        sg::solver::ConjugateGradients myCG(this->imax_final, this->epsilon_final);
        myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
      } else {
        sg::solver::ConjugateGradients myCG(this->imax, this->epsilon);
        myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
      }

      delete C;
    }

  }
}
