// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DMVECTORIZATIONPADDINGASSISTANT_HPP
#define DMVECTORIZATIONPADDINGASSISTANT_HPP

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>

#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

namespace parallel {

/**
 * This class implements some helper routines for padding
 * DataMatrix data containers such that their number of rows
 * is divisable by a number specified by a given vectorization type.
 */
class DMVectorizationPaddingAssistant {
 public:
  /**
   * determines the blocking vector-length for a given
   * vectorization mode, for double precision numbers.
   *
   * @param vecType selected vectorization mode
   *
   * @return the blocking/vector length
  */
  static size_t getVecWidth(VectorizationType& vecType);

  /**
   * determines the blocking vector-length for a given
   * vectorization mode, for single precision numbers.
   *
   * @param vecType selected vectorization mode
   *
   * @return the blocking/vector length
   */
  static size_t getVecWidthSP(VectorizationType& vecType);

  /**
   * Pads a DataMatrix object
   *
   * @param dataset dataset that should be padded
   * @param vecType vectorization which is used later on
   *
   * @return number of rows in the padded DataMatrix object
   */
  static size_t padDataset(SGPP::base::DataMatrix& dataset,
                           VectorizationType& vecType);

  /**
   * Pads a DataMatrixSP object
   *
   * @param dataset dataset that should be padded
   * @param vecType vectorization which is used later on
   *
   * @return number of rows in the padded DataMatrixSP object
   */
  static size_t padDataset(SGPP::base::DataMatrixSP& dataset,
                           VectorizationType vecType);
};

}

}

#endif