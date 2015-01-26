/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#include <sgpp/parallel/datadriven/algorithm/AlgorithmAdaBoostSPVectorizedIdentity.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMWeightMatrixSPVectorizedIdentity.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/solver/sle/ConjugateGradientsSP.hpp>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/PrecisionConverter.hpp>

namespace sg {
  namespace parallel {

    AlgorithmAdaBoostSPVectorizedIdentity::AlgorithmAdaBoostSPVectorizedIdentity(sg::base::Grid& SparseGrid, size_t gridType, int gridLevel, sg::base::DataMatrix& trainData, sg::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, size_t IMAX_final, double eps_final, double firstLabel, double secondLabel, double threshold, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, double percentOfAda, VectorizationType vecMode, size_t mode) : AlgorithmAdaBoostBase(SparseGrid, gridType, static_cast<sg::base::HashGenerator::level_t>(gridLevel), trainData, trainDataClass, NUM, lambda, IMAX, eps, IMAX_final, eps_final, firstLabel, secondLabel, threshold, maxLambda, minLambda, searchNum, refine, refineMode, refineNum, numberOfAda, percentOfAda, mode) {
      if (vecMode != X86SIMD && vecMode != OpenCL && vecMode != ArBB && vecMode != Hybrid_X86SIMD_OpenCL) {
        throw new sg::base::operation_exception("AlgorithmAdaBoostVectorizedIdentity : Only X86SIMD or OCL or ArBB or HYBRID_X86SIMD_OCL are supported vector extensions!");
      }

      this->vecMode = vecMode;
    }

    AlgorithmAdaBoostSPVectorizedIdentity::~AlgorithmAdaBoostSPVectorizedIdentity() {

    }

    void AlgorithmAdaBoostSPVectorizedIdentity::alphaSolver(double& lambda, sg::base::DataVector& weight, sg::base::DataVector& alpha, bool final) {
      sg::base::DataVectorSP classesSP(this->classes->getSize());
      sg::base::DataVectorSP alphaSP(alpha.getSize());
      sg::base::DataVectorSP weightSP(weight.getSize());
      sg::base::DataMatrixSP dataSP(this->data->getNrows(), this->data->getNcols());
      sg::base::DataVectorSP rhsSP(alpha.getSize());

      sg::base::PrecisionConverter::convertDataVectorToDataVectorSP(*this->classes, classesSP);
      sg::base::PrecisionConverter::convertDataVectorToDataVectorSP(alpha, alphaSP);
      sg::base::PrecisionConverter::convertDataVectorToDataVectorSP(weight, weightSP);
      sg::base::PrecisionConverter::convertDataMatrixToDataMatrixSP(*this->data, dataSP);

      DMWeightMatrixSPVectorizedIdentity WMatrix(*this->grid, dataSP, static_cast<float>(lambda), weightSP, this->vecMode);
      WMatrix.generateb(classesSP, rhsSP);

      if (final) {
        sg::solver::ConjugateGradientsSP myCG(this->imax_final, static_cast<float>(this->epsilon_final));
        myCG.solve(WMatrix, alphaSP, rhsSP, false, false, -1.0);
      } else {
        sg::solver::ConjugateGradientsSP myCG(this->imax, static_cast<float>(this->epsilon));
        myCG.solve(WMatrix, alphaSP, rhsSP, false, false, -1.0);
      }

      sg::base::PrecisionConverter::convertDataVectorSPToDataVector(alphaSP, alpha);
    }

  }

}
