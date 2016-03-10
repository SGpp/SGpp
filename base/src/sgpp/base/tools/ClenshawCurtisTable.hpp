// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CLENSHAWCURTISTABLE_HPP
#define CLENSHAWCURTISTABLE_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>

#include <cmath>
#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Lookup table for 1D Clenshaw-Curtis points.
 * This class precomputes the first \c maxLevel levels of a 1D Clenshaw-Curtis
 * grid to increase performance of Clenshaw-Curtis grids.
 */
class ClenshawCurtisTable {
 public:
  typedef HashGridIndex::level_type level_type;
  typedef HashGridIndex::index_type index_type;

  /// default number of intervals
  static const level_type DEFAULT_MAX_LEVEL = 16;

  static ClenshawCurtisTable& getInstance();

  /**
   * @param l       level of the grid point
   * @param i       index of the grid point (can be even)
   */
  inline double getPoint(level_type l, index_type i) const {
    if (l <= maxLevel) {
      return table.get((static_cast<index_type>(1) << l) +
                       static_cast<index_type>(l) + i - 1);
    } else {
      return calculatePoint(static_cast<index_type>(1) << l, i);
    }
  }

  /**
   * @param l       level of the grid point
   * @param i       index of the grid point (can be even)
   * @param hInv    2^l
   */
  inline double getPoint(level_type l, index_type i,
                          index_type hInv) const {
    if (l <= maxLevel) {
      return table.get(hInv +
                       static_cast<index_type>(l) + i - 1);
    } else {
      return calculatePoint(hInv, i);
    }
  }

 protected:
  /// lookup table
  DataVector table;
  /// maximal level
  level_type maxLevel;

  /**
   * @param hInv    2^l
   * @param i       index of the grid point (can be even)
   */
  inline double calculatePoint(index_type hInv, index_type i) const {
    return calculatePoint(
             1.0 / static_cast<double>(hInv), i);
  }

  /**
   * @param h       step width of the grid point (2^(-l))
   * @param i       index of the grid point (can be even)
   */
  inline double calculatePoint(double h, index_type i) const {
    return (std::cos(
              M_PI * (1.0 - static_cast<double>(i) * h)) + 1.0) / 2.0;
  }

 private:
  /**
   * Constructor creating the lookup table.
   *
   * @param maxLevel    level up to which grid points should be pre-computed
   */
  explicit ClenshawCurtisTable(level_type maxLevel = DEFAULT_MAX_LEVEL);

  /**
   * Deleted copy constructor.
   */
  ClenshawCurtisTable(const ClenshawCurtisTable&) = delete;

  /**
   * Deleted assignment operator.
   */
  void operator=(const ClenshawCurtisTable&) = delete;
};

}  // namespace base
}  // namespace sgpp

#endif /* CLENSHAWCURTISTABLE_HPP */
