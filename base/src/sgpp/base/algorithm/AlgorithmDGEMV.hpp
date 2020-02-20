// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ALGORTIHMDGEMV_HPP
#define ALGORTIHMDGEMV_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>
#include <utility>
#include <iostream>


namespace sgpp {
namespace base {

/**
 * Basic multiplaction with B and B^T on grids with no boundaries.
 * If there are @f$N@f$ basis functions @f$\varphi(\vec{x})@f$ and @f$m@f$ data points, then B is a (Nxm) matrix, with
 * @f[ (B)_{i,j} = \varphi_i(x_j). @f]
 * (The common known name for this operation is the BLAS routine DGEMV.)
 *
 */
template<class BASIS>
class AlgorithmDGEMV {
 public:
  /**
   * Performs the DGEMV Operation on the grid
   *
   * This operation can be executed in parallel by setting the USEOMP define
   *
   * @param storage GridStorage object that contains the grid's points information
   * @param basis a reference to a class that implements a specific basis
   * @param source the coefficients of the grid points
   * @param x the d-dimensional vector with data points (row-wise)
   * @param result the result vector of the matrix vector multiplication
   */
  void mult_transposed(GridStorage& storage, BASIS& basis,
                       const DataVector& source, DataMatrix& x, DataVector& result) {
    typedef std::vector<std::pair<size_t, double> > IndexValVector;

    result.setAll(0.0);

    #pragma omp parallel
    {
      size_t source_size = source.getSize();
      DataVector privateResult(result);
      DataVector line(x.getNcols());
      IndexValVector vec;
      GetAffectedBasisFunctions<BASIS> ga(storage);

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
  // implementation requires OpenMP 4.0 support
  //        void mult_transposed(GridStorage& storage, BASIS& basis,
  //        const DataVector& source, DataMatrix& x, DataVector& result) {
  //          typedef std::vector<std::pair<size_t, double> > IndexValVector;
  //
  //          result.setAll(0.0);
  //          #pragma omp declare reduction(accumulate :
  // ... sgpp::base::DataVector : omp_out.add(omp_in))
  // ... initializer ( omp_priv = DataVector(omp_orig.getSize(), 0))
  //
  //
  //          #pragma omp parallel
  //          {
  //            size_t source_size = source.getSize();
  //            DataVector line(x.getNcols());
  //            IndexValVector vec;
  //            GetAffectedBasisFunctions<BASIS> ga(storage);
  //
  //
  //            #pragma omp for reduction(accumulate:result) schedule(static)
  //
  //            for (size_t i = 0; i < source_size; i++) {
  //              vec.clear();
  //
  //              x.getRow(i, line);
  //
  //              ga(basis, line, vec);
  //
  //              for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
  //                result[iter->first] += iter->second * source[i];
  //              }
  //            }
  //          }
  //        }


  /**
   * Performs the DGEMV Operation on the grid having a transposed matrix
   *
   * This operation can be executed in parallel by setting the USEOMP define
   *
   *
   * @param storage GridStorage object that contains the grid's points information
   * @param basis a reference to a class that implements a specific basis
   * @param source the coefficients of the grid points
   * @param x the d-dimensional vector with data points (row-wise)
   * @param result the result vector of the matrix vector multiplication
   */
  void mult(GridStorage& storage, BASIS& basis, const DataVector& source,
            DataMatrix& x, DataVector& result) {
    typedef std::vector<std::pair<size_t, double> > IndexValVector;

    result.setAll(0.0);

    #pragma omp parallel
    {
      size_t result_size = result.getSize();

      DataVector line(x.getNcols());
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

}  // namespace base
}  // namespace sgpp

#endif /* ALGORTIHMDGEMV_HPP */
