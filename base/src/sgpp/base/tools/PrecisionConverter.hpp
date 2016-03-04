// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PRECISIONCONVERTER_HPP
#define PRECISIONCONVERTER_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Copies data from DataVector objects into DataVectorSP objects and vice versa.
 * Copies data from DataMatrix objects into DataMatrixSP objects and vice versa.
 *
 */
class PrecisionConverter {
 public:
  /**
   * Copies data from DataVector object into DataVectorSP object
   *
   * @param src source DataVector, double precision
   * @param dest destination DataVectorSP, single precision
   */
  static void convertDataVectorToDataVectorSP(const sgpp::base::DataVector& src,
      sgpp::base::DataVectorSP& dest);

  /**
   * Copies data from DataVectorSP object into DataVector object
   *
   * @param src source DataVectorSP, single precision
   * @param dest destination DataVector, double precision
   */
  static void convertDataVectorSPToDataVector(
    const sgpp::base::DataVectorSP& src,
    sgpp::base::DataVector& dest);

  /**
   * Copies data from DataMatrix object into DataMatrixSP object
   *
   * @param src source DataMatrix, double precision
   * @param dest destination DataMatrixSP, single precision
   */
  static void convertDataMatrixToDataMatrixSP(const sgpp::base::DataMatrix& src,
      sgpp::base::DataMatrixSP& dest);

  /**
   * Copies data from DataMatrixSP object into DataMatrix object
   *
   * @param src source DataMatrixSP, single precision
   * @param dest destination DataMatrix, double precision
   */
  static void convertDataMatrixSPToDataMatrix(
    const sgpp::base::DataMatrixSP& src,
    sgpp::base::DataMatrix& dest);
};

}  // namespace base
}  // namespace sgpp

#endif /* PRECISIONCONVERTER_HPP */
