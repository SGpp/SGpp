/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/operation/OperationEllipticPDESolverSystemDirichlet.hpp"
#include "base/exception/algorithm_exception.hpp"

#include <sstream>

namespace sg {
  namespace pde {

    OperationEllipticPDESolverSystemDirichlet::OperationEllipticPDESolverSystemDirichlet(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs) : OperationEllipticPDESolverSystem(SparseGrid, rhs) {
      this->BoundaryUpdate = new sg::base::DirichletUpdateVector(SparseGrid.getStorage());
      this->GridConverter = new sg::base::DirichletGridConverter();

      this->GridConverter->buildInnerGridWithCoefs(*this->BoundGrid, *this->rhs, &(this->InnerGrid), &(this->rhs_inner));

      this->numGridpointsInner = this->InnerGrid->getSize();

      this->alpha_inner = NULL;
    }

    OperationEllipticPDESolverSystemDirichlet::~OperationEllipticPDESolverSystemDirichlet() {
      delete this->alpha_inner;
      delete this->rhs_inner;
      delete this->InnerGrid;
      delete this->BoundaryUpdate;
      delete this->GridConverter;
    }

    void OperationEllipticPDESolverSystemDirichlet::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      applyLOperatorInner(alpha, result);
    }

    sg::base::DataVector* OperationEllipticPDESolverSystemDirichlet::generateRHS() {
      if (this->InnerGrid != NULL) {
        sg::base::DataVector alpha_tmp_complete(*(this->rhs));
        sg::base::DataVector rhs_tmp_complete(*(this->rhs));

        this->BoundaryUpdate->setInnerPointsToZero(alpha_tmp_complete);
        applyLOperatorComplete(alpha_tmp_complete, rhs_tmp_complete);

        this->GridConverter->calcInnerCoefs(rhs_tmp_complete, *(this->rhs_inner));
        this->rhs_inner->mult(-1.0);
      } else {
        throw new sg::base::algorithm_exception("OperationEllipticPDESolverSystemDirichlet::generateRHS : No inner grid exists!");
      }

      return this->rhs_inner;
    }

    sg::base::DataVector* OperationEllipticPDESolverSystemDirichlet::getGridCoefficientsForCG() {
      if (this->InnerGrid != NULL) {
        if (this->alpha_inner != NULL) {
          delete this->alpha_inner;
        }

        this->alpha_inner = new sg::base::DataVector(this->InnerGrid->getSize());
        this->alpha_inner->setAll(0.0);
      } else {
        throw new sg::base::algorithm_exception("OperationEllipticPDESolverSystemDirichlet::getGridCoefficientsForCG : No inner grid exists!");
      }

      return this->alpha_inner;
    }

    void OperationEllipticPDESolverSystemDirichlet::getSolutionBoundGrid(sg::base::DataVector& Solution, sg::base::DataVector& SolutionInner) {
      Solution = *(this->rhs);
      this->GridConverter->updateBoundaryCoefs(Solution, SolutionInner);
    }

    size_t OperationEllipticPDESolverSystemDirichlet::getMatrix(std::string& mtxString, bool complete) {
      size_t vector_size;

      if (complete == true) {
        vector_size = this->numGridpointsComplete;
      } else {
        vector_size = this->numGridpointsInner;
      }

      sg::base::DataVector alpha(vector_size);
      sg::base::DataVector result(vector_size);
      std::stringstream mtxStream;
      std::stringstream mtxHeaderStream;
      size_t nonZeros = 0;

      mtxHeaderStream.clear();
      mtxStream.clear();

      // loop over the system matrices columns
      for (size_t i = 0; i < vector_size; i++) {
        alpha.setAll(0.0);
        result.setAll(0.0);
        alpha.set(i, 1.0);

        // calculate column via Up/Down
        if (complete == true) {
          applyLOperatorComplete(alpha, result);
        } else {
          applyLOperatorInner(alpha, result);
        }

        // serialize result into mtxStream
        for (size_t j = 0; j < vector_size; j++) {
          if (result[j] != 0.0) {
            mtxStream << (j + 1) << " " << (i + 1) << " " << std::scientific << result[j] << std::endl;
            nonZeros++;
          }
        }
      }

      // Generate Header Line
      mtxHeaderStream << "%%MatrixMarket matrix coordinate real general" << std::endl << vector_size << " " << vector_size << " " << nonZeros << std::endl;

      mtxString = mtxHeaderStream.str() + mtxStream.str();

      return nonZeros;
    }

    void OperationEllipticPDESolverSystemDirichlet::getMatrixDiagonal(std::string& mtxString, bool complete) {
      size_t vector_size;

      if (complete == true) {
        vector_size = this->numGridpointsComplete;
      } else {
        vector_size = this->numGridpointsInner;
      }

      sg::base::DataVector alpha(vector_size);
      sg::base::DataVector result(vector_size);
      std::stringstream mtxStream;
      std::stringstream mtxHeaderStream;

      mtxHeaderStream.clear();
      mtxStream.clear();

      // loop over the system matrices columns
      for (size_t i = 0; i < vector_size; i++) {
        alpha.setAll(0.0);
        result.setAll(0.0);
        alpha.set(i, 1.0);

        // calculate column via Up/Down
        if (complete == true) {
          applyLOperatorComplete(alpha, result);
        } else {
          applyLOperatorInner(alpha, result);
        }

        // serialize result into mtxStream
        mtxStream << (i + 1) << " " << (i + 1) << " " << std::scientific << result[i] << std::endl;
      }

      // Generate Header Line
      mtxHeaderStream << "%%MatrixMarket matrix coordinate real general" << std::endl << vector_size << " " << vector_size << " " << vector_size << std::endl;

      mtxString = mtxHeaderStream.str() + mtxStream.str();
    }

    void OperationEllipticPDESolverSystemDirichlet::getMatrixDiagonalRowSum(std::string& mtxString, bool complete) {
      size_t vector_size;

      if (complete == true) {
        vector_size = this->numGridpointsComplete;
      } else {
        vector_size = this->numGridpointsInner;
      }

      sg::base::DataVector alpha(vector_size);
      sg::base::DataVector result(vector_size);
      sg::base::DataVector sum(vector_size);
      std::stringstream mtxStream;
      std::stringstream mtxHeaderStream;

      mtxHeaderStream.clear();
      mtxStream.clear();
      sum.setAll(0.0);

      // loop over the system matrices columns
      for (size_t i = 0; i < vector_size; i++) {
        alpha.setAll(0.0);
        result.setAll(0.0);
        alpha.set(i, 1.0);

        // calculate column via Up/Down
        if (complete == true) {
          applyLOperatorComplete(alpha, result);
        } else {
          applyLOperatorInner(alpha, result);
        }

        sum.add(result);
      }

      // generate Diagonalmatrix
      for (size_t i = 0; i < vector_size; i++) {
        // serialize result into mtxStream
        mtxStream << (i + 1) << " " << (i + 1) << " " << std::scientific << sum[i] << std::endl;
      }

      // Generate Header Line
      mtxHeaderStream << "%%MatrixMarket matrix coordinate real general" << std::endl << vector_size << " " << vector_size << " " << vector_size << std::endl;

      mtxString = mtxHeaderStream.str() + mtxStream.str();
    }

    size_t OperationEllipticPDESolverSystemDirichlet::getCompleteMatrix(std::string& mtxString) {
      return getMatrix(mtxString, true);
    }

    size_t OperationEllipticPDESolverSystemDirichlet::getInnerMatrix(std::string& mtxString) {
      return getMatrix(mtxString, false);
    }

    void OperationEllipticPDESolverSystemDirichlet::getInnerMatrixDiagonal(std::string& mtxString) {
      getMatrixDiagonal(mtxString, false);
    }

    void OperationEllipticPDESolverSystemDirichlet::getInnerMatrixDiagonalRowSum(std::string& mtxString) {
      getMatrixDiagonalRowSum(mtxString, false);
    }

  }
}
