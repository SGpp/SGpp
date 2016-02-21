// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ALGORITHMEVALUATION_HPP
#define ALGORITHMEVALUATION_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBoundaryBasis.hpp>

#include <sgpp/globaldef.hpp>

#include <utility>


namespace SGPP {
namespace base {

/**
 * Basic algorithm for getting all affected basis functions.
 * This implicitly assumes a tensor-product approach and local support.
 * No grid points on the border are supported.
 *
 * The main idea behind this algorithm is to spend as few function evaluations as possible.
 * Assume a regular sparse grid level 3 in two dimensions with the sparse grid basis
 * \f$\Phi:=\{\phi_i(x), i=1,\ldots,N\}\f$. Then the tableau of subspaces looks
 * as follows:
 *   \image html GetAffectedBasisFunctions_subspaces.png "Tableau of subspaces for a regular sparse grid level 3"
 * You could evaluate the function \f$ f_N(x) = \sum_{i=1}^N \alpha_i \phi_i(x)\f$ for all basis
 * functions \f$\phi_i(x)\f$, multiply them with the surplus and add them up.
 * In \f$d\f$ dimensions this would lead to \f$N\f$ evaluations of \f$d\f$ one-dimensional basis
 * functions each.
 *
 * A better way is to (recursively) look at each subspace, as only one basis function
 * per subspace can be non-zero (partially disjunct supports):
 *   \image html GetAffectedBasisFunctions_subspaces_affectedBasisFunctions.png "Traversal of subspaces for evaluation"
 * This can be done recursively in both the dimension and the level. In each subspace
 * the basis function concerned can be identified via a few index calculations and
 * evaluated at the given point in the domain.
 *
 * Even better would be to save further function evaluations and to reuse intermediate values obtained by
 * the evaluation of one-dimensional basis functions, see the following figure.
 *   \image html GetAffectedBasisFunctions_subspaces_affectedBasisFunctions_recursive.png "Minimize the number of evaluations" width=10cm
 * Descending recursively in the d-th dimension, one can propagate the value of the intermediate function
 * evaluation for the first d-1 dimensions that have already been looked at.
 *
 */
template<class BASIS>
class AlgorithmEvaluation {
 public:
  explicit AlgorithmEvaluation(GridStorage& storage) :
    storage(storage) {
  }

  ~AlgorithmEvaluation() {
  }

  /**
   * Returns evaluations of all basis functions that are non-zero at a given evaluation point.
   * For a given evaluation point \f$x\f$, it stores tuples (std::pair) of
   * \f$(i,\phi_i(x))\f$ in the result vector for all basis functions that are non-zero.
   * If one wants to evaluate \f$f_N(x)\f$, one only has to compute
   * \f[ \sum_{r\in\mathbf{result}} \alpha[r\rightarrow\mathbf{first}] \cdot r\rightarrow\mathbf{second}. \f]
   *
   * @param basis a sparse grid basis
   * @param point evaluation point within the domain
   * @param alpha the sparse grid's coefficients
   *
   * @result result result of the function evaluation
   */
  float_t operator()(BASIS& basis, const DataVector& point,
                     const DataVector& alpha) {
    GridStorage::grid_iterator working(storage);

    // typedef GridStorage::index_type::level_type level_type;
    typedef GridStorage::index_type::index_type index_type;

    size_t bits = sizeof(index_type) *
                  8;  // how many levels can we store in a index_type?

    size_t dim = storage.getDimension();

    // Check for bounding box
    BoundingBox* bb = storage.getBoundingBox();
    DataVector newPoint(point);

    if ( bb != NULL ) {
      for (size_t d = 0; d < dim; ++d) {
        DimensionBoundary dimbb = bb->getBoundary(d);
        // std::cout << "Dimension: " << d << " (left: " <<
        // dimbb.leftBoundary << ", right: " << dimbb.rightBoundary << ")" <<
        // std::endl;

        if (dimbb.leftBoundary == 0.0 && dimbb.rightBoundary == 1.0) {
          continue;
        }

        if (!(dimbb.leftBoundary <= newPoint[d]
              && newPoint[d] <= dimbb.rightBoundary) ) {
          // std::cout << "Out of bounds: " << point[d] << std::endl;
          return 0.0;
        }

        // std::cout << "Old: " << point[d] << std::endl;
        newPoint[d] = (newPoint[d] - dimbb.leftBoundary) /
                      (dimbb.rightBoundary - dimbb.leftBoundary);
        // std::cout << "New: " << point[d] << std::endl;
      }
    }

    index_type* source = new index_type[dim];

    for (size_t d = 0; d < dim; ++d) {
      // This does not really work on grids with borders.
      float_t temp = floor(newPoint[d] *
                           static_cast<float_t>(1 << (bits - 2))) * 2;

      if (newPoint[d] == 1.0) {
        source[d] = static_cast<index_type> (temp - 1);
      } else {
        source[d] = static_cast<index_type> (temp + 1);
      }
    }

    float_t result = 0.0;

    rec(basis, point, 0, 1.0, working, source, alpha, result);

    delete[] source;

    return result;
  }

 protected:
  GridStorage& storage;

  /**
   * Recursive traversal of the "tree" of basis functions for evaluation, used in operator().
   * For a given evaluation point \f$x\f$, it stores tuples (std::pair) of
   * \f$(i,\phi_i(x))\f$ in the result vector for all basis functions that are non-zero.
   *
   * @param basis a sparse grid basis
   * @param point evaluation point within the domain
   * @param current_dim the dimension currently looked at (recursion parameter)
   * @param value the value of the evaluation of the current basis function up to (excluding) dimension current_dim (product of the evaluations of the one-dimensional ones)
   * @param working iterator working on the GridStorage of the basis
   * @param source array of indices for each dimension (identifying the indices of the current grid point)
   * @param alpha the spars grid's ansatzfunctions coefficients
   * @param result reference to a float_t into which the result should be stored
   */
  void rec(BASIS& basis, const DataVector& point, size_t current_dim,
           float_t value, GridStorage::grid_iterator& working,
           GridStorage::index_type::index_type* source, const DataVector& alpha,
           float_t& result) {
    typedef GridStorage::index_type::level_type level_type;
    typedef GridStorage::index_type::index_type index_type;

    const unsigned int BITS_IN_BYTE = 8;
    // maximum possible level for the index type
    const level_type max_level = static_cast<level_type> (
                                   sizeof(index_type) * BITS_IN_BYTE - 1);
    index_type src_index = source[current_dim];

    level_type work_level = 1;

    while (true) {
      size_t seq = working.seq();

      if (storage.end(seq)) {
        break;
      } else {
        index_type work_index;
        level_type temp;

        working.get(current_dim, temp, work_index);

        float_t new_value = basis.eval(work_level, work_index,
                                       point[current_dim]);
        new_value *= value;

        if (current_dim == storage.getDimension() - 1) {
          result += (alpha[seq] * new_value);
        } else {
          rec(basis, point, current_dim + 1, new_value,
              working, source, alpha, result);
        }
      }

      if (working.hint()) {
        break;
      }

      // this decides in which direction we should descend by evaluating
      // the corresponding bit
      // the bits are coded from left to right starting with level 1
      // being in position max_level
      bool right = (src_index & (1 << (max_level - work_level))) > 0;
      ++work_level;

      if (right) {
        working.rightChild(current_dim);
      } else {
        working.leftChild(current_dim);
      }
    }

    working.resetToLevelOne(current_dim);
  }
};

}  // namespace base
}  // namespace SGPP

#endif /* ALGORITHMEVALUATION_HPP */
