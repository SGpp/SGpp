/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/tools/PrecisionConverter.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void PrecisionConverter::convertDataVectorToDataVectorSP(const SGPP::base::DataVector& src, SGPP::base::DataVectorSP& dest) {
      if (src.getSize() != dest.getSize()) {
        throw new SGPP::base::data_exception("PrecisionConverter::convertDataVectorToDataVectorSP : vector sizes don't match!");
      } else {
        for (size_t i = 0; i < src.getSize(); i++) {
          dest.set(i, static_cast<float>(src.get(i)));
        }
      }
    }

    void PrecisionConverter::convertDataVectorSPToDataVector(const SGPP::base::DataVectorSP& src, SGPP::base::DataVector& dest) {
      if (src.getSize() != dest.getSize()) {
        throw new SGPP::base::data_exception("PrecisionConverter::convertDataVectorSPToDataVector : vector sizes don't match!");
      } else {
        for (size_t i = 0; i < src.getSize(); i++) {
          dest.set(i, static_cast<double>(src.get(i)));
        }
      }
    }

    void PrecisionConverter::convertDataMatrixToDataMatrixSP(const SGPP::base::DataMatrix& src, SGPP::base::DataMatrixSP& dest) {
      if (src.getNcols() != dest.getNcols() || src.getNrows() != dest.getNrows()) {
        throw new SGPP::base::data_exception("PrecisionConverter::convertDataMatrixToDataMatrixSP : matrix sizes don't match!");
      } else {
        for (size_t i = 0; i < src.getNrows(); i++) {
          for (size_t j = 0; j < src.getNcols(); j++) {
            dest.set(i, j, static_cast<float>(src.get(i, j)));
          }
        }
      }
    }

    void PrecisionConverter::convertDataMatrixSPToDataMatrix(const SGPP::base::DataMatrixSP& src, SGPP::base::DataMatrix& dest) {
      if (src.getNcols() != dest.getNcols() || src.getNrows() != dest.getNrows()) {
        throw new SGPP::base::data_exception("PrecisionConverter::convertDataMatrixSPToDataMatrix : matrix sizes don't match!");
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
