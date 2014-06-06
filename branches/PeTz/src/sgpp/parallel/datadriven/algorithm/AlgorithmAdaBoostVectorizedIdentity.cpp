/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Zhongwen Song (songz@in.tum.de)
// @author Benjamin (pehersto@in.tum.de)
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#include "parallel/datadriven/algorithm/AlgorithmAdaBoostVectorizedIdentity.hpp"
#include "parallel/datadriven/algorithm/DMWeightMatrixVectorizedIdentity.hpp"
#include "base/exception/operation_exception.hpp"
#include "base/operation/BaseOpFactory.hpp"

namespace sg {
  namespace parallel {

    AlgorithmAdaBoostVectorizedIdentity::AlgorithmAdaBoostVectorizedIdentity(sg::base::Grid& SparseGrid, size_t gridType, int gridLevel, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, size_t IMAX_final, double eps_final, double firstLabel, double secondLabel, double threshold, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, double percentOfAda, VectorizationType vecMode, size_t mode)
      : AlgorithmAdaBoostBase(SparseGrid, gridType, static_cast<sg::base::HashGenerator::level_t>(gridLevel), trainData, trainDataClass, NUM, lambda, IMAX, eps, IMAX_final, eps_final, firstLabel, secondLabel, threshold, maxLambda, minLambda, searchNum, refine, refineMode, refineNum, numberOfAda, percentOfAda, mode) {
      if (vecMode != X86SIMD && vecMode != OpenCL && vecMode != ArBB && vecMode != Hybrid_X86SIMD_OpenCL) {
        throw new sg::base::operation_exception("AlgorithmAdaBoostVectorizedIdentity : Only X86SIMD or OCL or ArBB or HYBRID_X86SIMD_OCL are supported vector extensions!");
      }

      this->vecMode = vecMode;
    }

    AlgorithmAdaBoostVectorizedIdentity::~AlgorithmAdaBoostVectorizedIdentity() {

    }

    void AlgorithmAdaBoostVectorizedIdentity::alphaSolver(double& lambda, sg::base::DataVector& weight, sg::base::DataVector& alpha, bool final) {
      DMWeightMatrixVectorizedIdentity WMatrix(*this->grid, *this->data, lambda, weight, this->vecMode);
      sg::base::DataVector rhs(alpha.getSize());
      WMatrix.generateb(*this->classes, rhs);

      if (final) {
        sg::solver::ConjugateGradients myCG(this->imax_final, this->epsilon_final);
        myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
      } else {
        sg::solver::ConjugateGradients myCG(this->imax, this->epsilon);
        myCG.solve(WMatrix, alpha, rhs, false, false, -1.0);
      }
    }

  }

}
