/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef PRECISIONCONVERTER_HPP
#define PRECISIONCONVERTER_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/datatypes/DataVectorSP.hpp"
#include "base/datatypes/DataMatrixSP.hpp"


namespace sg {
  namespace base {

    /**
     * Copies data from DataVector objects into DataVectorSP objects and vice versa.
     * Copies data from DataMatrix objects into DataMatrixSP objects and vice versa.
     *
     * @version $HEAD$
     */
    class PrecisionConverter {
      public:
        /**
         * Copies data from DataVector object into DataVectorSP object
         *
         * @param src source DataVector, double precision
         * @param dest destination DataVectorSP, single precision
         */
        static void convertDataVectorToDataVectorSP(const sg::base::DataVector& src, sg::base::DataVectorSP& dest);

        /**
         * Copies data from DataVectorSP object into DataVector object
         *
         * @param src source DataVectorSP, single precision
         * @param dest destination DataVector, double precision
         */
        static void convertDataVectorSPToDataVector(const sg::base::DataVectorSP& src, sg::base::DataVector& dest);

        /**
         * Copies data from DataMatrix object into DataMatrixSP object
         *
         * @param src source DataMatrix, double precision
         * @param dest destination DataMatrixSP, single precision
         */
        static void convertDataMatrixToDataMatrixSP(const sg::base::DataMatrix& src, sg::base::DataMatrixSP& dest);

        /**
         * Copies data from DataMatrixSP object into DataMatrix object
         *
         * @param src source DataMatrixSP, single precision
         * @param dest destination DataMatrix, double precision
         */
        static void convertDataMatrixSPToDataMatrix(const sg::base::DataMatrixSP& src, sg::base::DataMatrix& dest);
    };

  }

}

#endif /* PRECISIONCONVERTER_HPP */
