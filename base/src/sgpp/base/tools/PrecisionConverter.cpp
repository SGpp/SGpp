/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/tools/PrecisionConverter.hpp>
#include <sgpp/base/exception/data_exception.hpp>

namespace sg {
  namespace base {

    void PrecisionConverter::convertDataVectorToDataVectorSP(const sg::base::DataVector& src, sg::base::DataVectorSP& dest) {
      if (src.getSize() != dest.getSize()) {
        throw new sg::base::data_exception("PrecisionConverter::convertDataVectorToDataVectorSP : vector sizes don't match!");
      } else {
        for (size_t i = 0; i < src.getSize(); i++) {
          dest.set(i, static_cast<float>(src.get(i)));
        }
      }
    }

    void PrecisionConverter::convertDataVectorSPToDataVector(const sg::base::DataVectorSP& src, sg::base::DataVector& dest) {
      if (src.getSize() != dest.getSize()) {
        throw new sg::base::data_exception("PrecisionConverter::convertDataVectorSPToDataVector : vector sizes don't match!");
      } else {
        for (size_t i = 0; i < src.getSize(); i++) {
          dest.set(i, static_cast<double>(src.get(i)));
        }
      }
    }

    void PrecisionConverter::convertDataMatrixToDataMatrixSP(const sg::base::DataMatrix& src, sg::base::DataMatrixSP& dest) {
      if (src.getNcols() != dest.getNcols() || src.getNrows() != dest.getNrows()) {
        throw new sg::base::data_exception("PrecisionConverter::convertDataMatrixToDataMatrixSP : matrix sizes don't match!");
      } else {
        for (size_t i = 0; i < src.getNrows(); i++) {
          for (size_t j = 0; j < src.getNcols(); j++) {
            dest.set(i, j, static_cast<float>(src.get(i, j)));
          }
        }
      }
    }

    void PrecisionConverter::convertDataMatrixSPToDataMatrix(const sg::base::DataMatrixSP& src, sg::base::DataMatrix& dest) {
      if (src.getNcols() != dest.getNcols() || src.getNrows() != dest.getNrows()) {
        throw new sg::base::data_exception("PrecisionConverter::convertDataMatrixSPToDataMatrix : matrix sizes don't match!");
      } else {
        for (size_t i = 0; i < src.getNrows(); i++) {
          for (size_t j = 0; j < src.getNcols(); j++) {
            dest.set(i, j, static_cast<double>(src.get(i, j)));
          }
        }
      }
    }

  }

}
