// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/datadriven/algorithm/AlgorithmAdaBoostSPVectorizedIdentity.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMWeightMatrixSPVectorizedIdentity.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/solver/sle/ConjugateGradientsSP.hpp>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/PrecisionConverter.hpp>

#include <sgpp/globaldef.hpp>

#if USE_DOUBLE_PRECISION==0


namespace SGPP {
  namespace parallel {

    AlgorithmAdaBoostSPVectorizedIdentity::AlgorithmAdaBoostSPVectorizedIdentity(SGPP::base::Grid& SparseGrid, size_t gridType, int gridLevel, SGPP::base::DataMatrix& trainData, SGPP::base::DataVector& trainDataClass, size_t NUM, double lambda, size_t IMAX, double eps, size_t IMAX_final, double eps_final, double firstLabel, double secondLabel, double threshold, double maxLambda, double minLambda, size_t searchNum, bool refine, size_t refineMode, size_t refineNum, size_t numberOfAda, double percentOfAda, VectorizationType vecMode, size_t mode) : AlgorithmAdaBoostBase(SparseGrid, gridType, static_cast<SGPP::base::HashGenerator::level_t>(gridLevel), trainData, trainDataClass, NUM, lambda, IMAX, eps, IMAX_final, eps_final, firstLabel, secondLabel, threshold, maxLambda, minLambda, searchNum, refine, refineMode, refineNum, numberOfAda, percentOfAda, mode) {
      if (vecMode != X86SIMD && vecMode != OpenCL && vecMode != ArBB && vecMode != Hybrid_X86SIMD_OpenCL) {
        throw new SGPP::base::operation_exception("AlgorithmAdaBoostVectorizedIdentity : Only X86SIMD or OCL or ArBB or HYBRID_X86SIMD_OCL are supported vector extensions!");
      }

      this->vecMode = vecMode;
    }

    AlgorithmAdaBoostSPVectorizedIdentity::~AlgorithmAdaBoostSPVectorizedIdentity() {

    }

    void AlgorithmAdaBoostSPVectorizedIdentity::alphaSolver(double& lambda, SGPP::base::DataVector& weight, SGPP::base::DataVector& alpha, bool final) {
      SGPP::base::DataVectorSP classesSP(this->classes->getSize());
      SGPP::base::DataVectorSP alphaSP(alpha.getSize());
      SGPP::base::DataVectorSP weightSP(weight.getSize());
      SGPP::base::DataMatrixSP dataSP(this->data->getNrows(), this->data->getNcols());
      SGPP::base::DataVectorSP rhsSP(alpha.getSize());

      SGPP::base::PrecisionConverter::convertDataVectorToDataVectorSP(*this->classes, classesSP);
      SGPP::base::PrecisionConverter::convertDataVectorToDataVectorSP(alpha, alphaSP);
      SGPP::base::PrecisionConverter::convertDataVectorToDataVectorSP(weight, weightSP);
      SGPP::base::PrecisionConverter::convertDataMatrixToDataMatrixSP(*this->data, dataSP);

      DMWeightMatrixSPVectorizedIdentity WMatrix(*this->grid, dataSP, static_cast<float>(lambda), weightSP, this->vecMode);
      WMatrix.generateb(classesSP, rhsSP);

      if (final) {
        SGPP::solver::ConjugateGradientsSP myCG(this->imax_final, static_cast<float>(this->epsilon_final));
        myCG.solve(WMatrix, alphaSP, rhsSP, false, false, -1.0);
      } else {
        SGPP::solver::ConjugateGradientsSP myCG(this->imax, static_cast<float>(this->epsilon));
        myCG.solve(WMatrix, alphaSP, rhsSP, false, false, -1.0);
      }

      SGPP::base::PrecisionConverter::convertDataVectorSPToDataVector(alphaSP, alpha);
    }

  }

}

#endif
