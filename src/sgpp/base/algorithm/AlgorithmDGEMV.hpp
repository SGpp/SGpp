/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ALGORTIHMDGEMV_HPP
#define ALGORTIHMDGEMV_HPP

#include "base/grid/GridStorage.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#include "base/algorithm/GetAffectedBasisFunctions.hpp"

#include <vector>
#include <utility>
#include <iostream>

namespace sg {
  namespace base {

    /**
     * Basic multiplaction with B and B^T on grids with no boundaries.
     * If there are @f$N@f$ basis functions @f$\varphi(\vec{x})@f$ and @f$m@f$ data points, then B is a (Nxm) matrix, with
     * @f[ (B)_{i,j} = \varphi_i(x_j). @f]
     * (The common known name for this operation is the BLAS routine DGEMV.)
     *
     * @todo (blank) check if it is possible to have some functor for the BASIS type
     */
    template<class BASIS>
    class AlgorithmDGEMV {
      public:

        /**
         * Performs the DGEMV Operation on the grid
         *
         * This operation can be executed in parallel by setting the USEOMP define
         *
         * @todo (heinecke, nice) add mathematical description
         *
         * @param storage GridStorage object that contains the grid's points information
         * @param basis a reference to a class that implements a specific basis
         * @param source the coefficients of the grid points
         * @param x the d-dimensional vector with data points (row-wise)
         * @param result the result vector of the matrix vector multiplication
         */
        void mult_transposed(GridStorage* storage, BASIS& basis, DataVector& source, DataMatrix& x, DataVector& result) {
          typedef std::vector<std::pair<size_t, double> > IndexValVector;

          result.setAll(0.0);

          #pragma omp parallel
          {
            size_t source_size = source.getSize();
            DataVector privateResult(result);
            std::vector<double> line;
            IndexValVector vec;
            GetAffectedBasisFunctions<BASIS> ga(storage);

            privateResult.setAll(0.0);

            #pragma omp for schedule(static)

            for (size_t i = 0; i < source_size; i++) {
              vec.clear();

              x.getRow(i, line);

              ga(basis, line, vec);

              for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
                privateResult[iter->first] += iter->second * source[i];
              }
            }

            #pragma omp critical
            {
              result.add(privateResult);
            }
          }
        }

        /**
         * Performs the DGEMV Operation on the grid having a transposed matrix
         *
         * This operation can be executed in parallel by setting the USEOMP define
         *
         * @todo (heinecke, nice) add mathematical description
         *
         * @param storage GridStorage object that contains the grid's points information
         * @param basis a reference to a class that implements a specific basis
         * @param source the coefficients of the grid points
         * @param x the d-dimensional vector with data points (row-wise)
         * @param result the result vector of the matrix vector multiplication
         */
        void mult(GridStorage* storage, BASIS& basis, DataVector& source, DataMatrix& x, DataVector& result) {
          typedef std::vector<std::pair<size_t, double> > IndexValVector;

          result.setAll(0.0);

          #pragma omp parallel
          {
            size_t result_size = result.getSize();

            std::vector<double> line;
            IndexValVector vec;

            GetAffectedBasisFunctions<BASIS> ga(storage);

            #pragma omp for schedule (static)

            for (size_t i = 0; i < result_size; i++) {
              vec.clear();

              x.getRow(i, line);

              ga(basis, line, vec);

              for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
                result[i] += iter->second * source[iter->first];
              }
            }
          }
        }
    };

  }
}

#endif /* ALGORTIHMDGEMV_HPP */
