/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef PRECISIONCONVERTER_HPP
#define PRECISIONCONVERTER_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVectorSP.hpp>
#include <sgpp/base/datatypes/DataMatrixSP.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
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
        static void convertDataVectorToDataVectorSP(const SGPP::base::DataVector& src, SGPP::base::DataVectorSP& dest);

        /**
         * Copies data from DataVectorSP object into DataVector object
         *
         * @param src source DataVectorSP, single precision
         * @param dest destination DataVector, double precision
         */
        static void convertDataVectorSPToDataVector(const SGPP::base::DataVectorSP& src, SGPP::base::DataVector& dest);

        /**
         * Copies data from DataMatrix object into DataMatrixSP object
         *
         * @param src source DataMatrix, double precision
         * @param dest destination DataMatrixSP, single precision
         */
        static void convertDataMatrixToDataMatrixSP(const SGPP::base::DataMatrix& src, SGPP::base::DataMatrixSP& dest);

        /**
         * Copies data from DataMatrixSP object into DataMatrix object
         *
         * @param src source DataMatrixSP, single precision
         * @param dest destination DataMatrix, double precision
         */
        static void convertDataMatrixSPToDataMatrix(const SGPP::base::DataMatrixSP& src, SGPP::base::DataMatrix& dest);
    };

  }

}

#endif /* PRECISIONCONVERTER_HPP */
