// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Potentially hetereogeneous basis on full grids, i.e., a <tt>dim</tt>-dimensional vector of
 * 1D bases of type sgpp::base::Basis. In every dimension, a different basis may be used
 * (or the same basis may be used for all dimensions).
 */
class HeterogeneousBasis {
 public:
  /**
   * Default constructor, corresponds to the zero-dimensional case.
   */
  HeterogeneousBasis() : bases1d() {
  }

  /**
   * Constructor.
   *
   * @param bases1d         vector of unique_ptr to 1D bases (do not destruct before this object)
   * @param isHierarchical  whether the basis is hierarchical or nodal
   */
  explicit HeterogeneousBasis(
      const std::vector<std::unique_ptr<base::Basis<level_t, index_t>>>& bases1d,
      bool isHierarchical = true) : bases1d(), isHierarchical_(isHierarchical) {
    for (const std::unique_ptr<base::Basis<level_t, index_t>>& basis1d : bases1d) {
      this->bases1d.push_back(basis1d.get());
    }
  }

  /**
   * Constructor.
   *
   * @param bases1d         vector of pointers to 1D bases (do not delete before this object)
   * @param isHierarchical  whether the basis is hierarchical or nodal
   */
  explicit HeterogeneousBasis(const std::vector<base::Basis<level_t, index_t>*>& bases1d,
      bool isHierarchical = true) : bases1d(bases1d), isHierarchical_(isHierarchical) {
  }

  /**
   * Constructor for the case in which the same basis should be used for all dimensions.
   *
   * @param dim             dimensionality
   * @param basis1d         1D basis (do not destruct before this object)
   * @param isHierarchical  whether the basis is hierarchical or nodal
   */
  HeterogeneousBasis(size_t dim, base::Basis<level_t, index_t>& basis1d,
      bool isHierarchical = true) : bases1d(dim, &basis1d), isHierarchical_(isHierarchical) {
  }

  /**
   * Equality operator. Two HeterogeneousBasis instances are equal if they contain the same
   * 1D basis pointers and if the values of \c isHierarchical are the same.
   *
   * @param other   other HeterogeneousBasis instance
   * @return whether both instances are equal
   */
  bool operator==(const HeterogeneousBasis& other) const {
    return (bases1d == other.bases1d) && (isHierarchical_ == other.isHierarchical_);
  }

  /**
   * Inequality operator. Two HeterogeneousBasis instances are equal if they contain the same
   * 1D basis pointers and if the values of \c isHierarchical are the same.
   *
   * @param other   other HeterogeneousBasis instance
   * @return whether both instances are inequal
   */
  bool operator!=(const HeterogeneousBasis& other) const {
    return !(*this == other);
  }

  /**
   * Hierarchize level/index pair. For every level \f$\ell \ge 0\f$ and index
   * \f$i = 0, \dotsc, 2^\ell\f$, there is a unique level \f$0 \le \tilde{\ell} \le \ell\f$ and
   * index \f$\tilde{i} = 0, \dotsc, 2^\tilde{\ell}\f$ such that
   *
   * - \f$i 2^{-\ell} =: x_{l,i} = x_{\tilde{l},\tilde{i}} := \tilde{i} 2^{-\tilde{\ell}}\f$ and
   * - \f$\tilde{i}\f$ is odd or \f$\tilde{l} = 0\f$.
   *
   * @param[in,out] level   \f$\ell\f$ before calling, \f$\tilde{\ell}\f$ after calling
   * @param[in,out] index   \f$i\f$ before calling, \f$\tilde{i}\f$ after calling
   */
  static void hierarchizeLevelIndex(level_t& level, index_t& index) {
    if (index == 0) {
      level = 0;
    } else if (level > 0) {
      while (index % 2 == 0) {
        level--;
        index /= 2;
      }
    }
  }

  /**
   * Hierarchize pair of level/index multi-indices. Applies 1D version to each dimension.
   *
   * @param[in,out] level   \f$\vec{\ell}\f$ before calling, \f$\tilde{\vec{\ell}}\f$ after calling
   * @param[in,out] index   \f$\vec{i}\f$ before calling, \f$\vec{\tilde{i}}\f$ after calling
   */
  static void hierarchizeLevelIndex(LevelVector& level, IndexVector& index) {
    for (size_t d = 0; d < level.size(); d++) {
      hierarchizeLevelIndex(level[d], index[d]);
    }
  }

  /**
   * Hierarchize vector of level/index multi-indices. Applies multi-index version to every entry.
   *
   * @param[in,out] levels      vector of \f$\vec{\ell}\f$ before calling,
   *                            vector of \f$\tilde{\vec{\ell}}\f$ after calling
   * @param[in,out] indices     vector of \f$\vec{i}\f$ before calling,
   *                            vector of \f$\vec{\tilde{i}}\f$ after calling
   */
  static void hierarchizeLevelIndex(std::vector<LevelVector>& levels,
      std::vector<IndexVector>& indices) {
    for (size_t i = 0; i < levels.size(); i++) {
      hierarchizeLevelIndex(levels[i], indices[i]);
    }
  }

  /**
   * Evaluate basis function of given level and index at given point as the tensor product of the
   * 1D bases. If the basis is hierarchical, \c hierarchizeLevelIndex will be called before
   * evaluating the 1D bases.
   *
   * @param level   level of the basis function
   * @param index   index of the basis function
   * @param point   point at which to evaluate the basis function
   * @return value of the basis function at the given point
   */
  double eval(const LevelVector& level, const IndexVector& index,
      const base::DataVector& point) const {
    double result = 1.0;
    const size_t dim = bases1d.size();

    if (isHierarchical_) {
      for (size_t d = 0; d < dim; d++) {
        level_t l = level[d];
        index_t i = index[d];
        hierarchizeLevelIndex(l, i);
        result *= bases1d[d]->eval(l, i, point[d]);
      }
    } else {
      for (size_t d = 0; d < dim; d++) {
        result *= bases1d[d]->eval(level[d], index[d], point[d]);
      }
    }

    return result;
  }

  /**
   * @return dimensionality
   */
  size_t getDimension() const {
    return bases1d.size();
  }

  /**
   * @return vector of pointers to 1D bases (do not delete before this object)
   */
  const std::vector<base::Basis<level_t, index_t>*>& getBases1d() const {
    return bases1d;
  }

  /**
   * @param bases1d   vector of pointers to 1D bases (do not delete before this object)
   */
  void setBases1D(const std::vector<base::Basis<level_t, index_t>*>& bases1d) {
    this->bases1d = bases1d;
  }

  /**
   * @return whether the basis is hierarchical or nodal
   */
  bool isHierarchical() const {
    return isHierarchical_;
  }

  /**
   * @param isHierarchical  whether the basis is hierarchical or nodal
   */
  void setIsHierarchical(bool isHierarchical) {
    isHierarchical_ = isHierarchical;
  }

 protected:
  /// vector of pointers to 1D bases (do not delete before this object)
  std::vector<base::Basis<level_t, index_t>*> bases1d;
  /// whether the basis is hierarchical or nodal
  bool isHierarchical_;
};

}  // namespace combigrid
}  // namespace sgpp
