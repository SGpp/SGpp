/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "datadriven/algorithm/DMSystemMatrixExplicit.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/operation/OperationMatrix.hpp"
#include <limits.h>

using namespace sg::base;

namespace sg {
  namespace datadriven {

    DMSystemMatrixExplicit::DMSystemMatrixExplicit( DataMatrix& trainData, double lambda)
      : DMSystemMatrixBase(trainData, lambda) {
      //this->dataset_->transpose();
    }

    DMSystemMatrixExplicit::~DMSystemMatrixExplicit()
    {}

    void DMSystemMatrixExplicit::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      //TODO: implement the method
    }

    void DMSystemMatrixExplicit::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
      //TODO: impplement the method
    }

    double* DMSystemMatrixExplicit::compileMatrix(sg::base::GridStorage* storage,unsigned int first_index,unsigned int last_index) {
      // Right now only support for modlin basis without vectorization is implemented
      //size_t source_size = this->dataset_->getNcols();
      size_t source_size = this->dataset_->getNrows();
      size_t dims = storage->dim();
      size_t storageSize = storage->size();
      last_index = (last_index == UINT_MAX) ? storageSize - 1 : last_index;
      sg::base::DataMatrix level(storage->size(), storage->dim());
      sg::base::DataMatrix index(storage->size(), storage->dim());
      storage->getLevelIndexArraysForEval(level, index);
      double* ptrData = this->dataset_->getPointer();
      double* ptrLevel = level.getPointer();
      double* ptrIndex = index.getPointer();
      double* ptrResult = new double[(last_index - first_index + 1)*source_size];
      int matrixIndex = 0;

      for (size_t i = 0; i < source_size; i++) {
        for (size_t j = first_index; j <= last_index; j++) {
          double curSupport = 1.0;

          for (size_t d = 0; d < dims; d++) {
            int data_index = (i * dims) + d; //(d*source_size)+i

            if (ptrLevel[(j * dims) + d] == 2.0) {
              // nothing to do (mult with 1)
            } else if (ptrIndex[(j * dims) + d] == 1.0) {
              //double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
              double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[data_index]));
              eval = 2.0 - eval;
              double localSupport = std::max<double>(eval, 0.0);
              curSupport *= localSupport;
            } else if (ptrIndex[(j * dims) + d] == (ptrLevel[(j * dims) + d] - 1.0)) {
              //double eval = ((ptrLevel[(j*dims)+d]) * (ptrData[(d*source_size)+i]));
              double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[data_index]));
              double index_calc = eval - (ptrIndex[(j * dims) + d]);
              double last = 1.0 + index_calc;
              double localSupport = std::max<double>(last, 0.0);
              curSupport *= localSupport;
            } else {
              double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[data_index]));
              double index_calc = eval - (ptrIndex[(j * dims) + d]);
              double abs = fabs(index_calc);
              double last = 1.0 - abs;
              double localSupport = std::max<double>(last, 0.0);
              curSupport *= localSupport;
            }
          }

          ptrResult[matrixIndex++] = curSupport;
        }
      }

      return ptrResult;
    }


  }
}
