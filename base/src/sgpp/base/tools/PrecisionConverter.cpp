// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/PrecisionConverter.hpp>
#include <sgpp/base/exception/data_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

void PrecisionConverter::convertDataVectorToDataVectorSP(
  const sgpp::base::DataVector& src, sgpp::base::DataVectorSP& dest) {
  if (src.getSize() != dest.getSize()) {
    throw sgpp::base::data_exception(
      "PrecisionConverter::convertDataVectorToDataVectorSP : "
      "vector sizes don't match!");
  } else {
    for (size_t i = 0; i < src.getSize(); i++) {
      dest.set(i, static_cast<float>(src.get(i)));
    }
  }
}

void PrecisionConverter::convertDataVectorSPToDataVector(
  const sgpp::base::DataVectorSP& src, sgpp::base::DataVector& dest) {
  if (src.getSize() != dest.getSize()) {
    throw sgpp::base::data_exception(
      "PrecisionConverter::convertDataVectorSPToDataVector : "
      "vector sizes don't match!");
  } else {
    for (size_t i = 0; i < src.getSize(); i++) {
      dest.set(i, static_cast<double>(src.get(i)));
    }
  }
}

void PrecisionConverter::convertDataMatrixToDataMatrixSP(
  const sgpp::base::DataMatrix& src, sgpp::base::DataMatrixSP& dest) {
  if (src.getNcols() != dest.getNcols() || src.getNrows() != dest.getNrows()) {
    throw sgpp::base::data_exception(
      "PrecisionConverter::convertDataMatrixToDataMatrixSP : "
      "matrix sizes don't match!");
  } else {
    for (size_t i = 0; i < src.getNrows(); i++) {
      for (size_t j = 0; j < src.getNcols(); j++) {
        dest.set(i, j, static_cast<float>(src.get(i, j)));
      }
    }
  }
}

void PrecisionConverter::convertDataMatrixSPToDataMatrix(
  const sgpp::base::DataMatrixSP& src, sgpp::base::DataMatrix& dest) {
  if (src.getNcols() != dest.getNcols() || src.getNrows() != dest.getNrows()) {
    throw sgpp::base::data_exception(
      "PrecisionConverter::convertDataMatrixSPToDataMatrix : "
      "matrix sizes don't match!");
  } else {
    for (size_t i = 0; i < src.getNrows(); i++) {
      for (size_t j = 0; j < src.getNcols(); j++) {
        dest.set(i, j, static_cast<double>(src.get(i, j)));
      }
    }
  }
}

}  // namespace base
}  // namespace sgpp
