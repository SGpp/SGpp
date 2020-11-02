// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/common/basis/Basis.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/system/SLE.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/basis/HeterogeneousBasis.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Operation of hierarchising values on a 1D pole of a full grid with general basis functions
 * by solving a linear system.
 */
class OperationPoleHierarchisationGeneral : public OperationPole {
 public:
  /**
   * Constructor.
   *
   * @param basis                 1D basis
   * @param isBasisHierarchical   whether the basis is hierarchical or nodal
   */
  explicit OperationPoleHierarchisationGeneral(base::Basis<level_t, index_t>& basis,
      bool isBasisHierarchical = true);

  /**
   * Virtual destructor.
   */
  ~OperationPoleHierarchisationGeneral() override;

  /**
   * Factory method to create a vector of unique_ptr to OperationPoleHierarchisationGeneral objects
   * from a HeterogeneousBasis.
   *
   * @param[in] basis       multivariate basis
   * @param[out] operation  vector of unique_ptr to OperationPoleHierarchisationGeneral objects
   *                        corresponding to the given basis
   */
  static void fromHeterogenerousBasis(const HeterogeneousBasis& basis,
      std::vector<std::unique_ptr<OperationPole>>& operation);

  /**
   * Factory method to create a vector of pointers to OperationPoleHierarchisationGeneral objects
   * from a HeterogeneousBasis.
   *
   * @param[in] basis       multivariate basis
   * @param[out] operation  vector of pointers to OperationPoleHierarchisationGeneral objects
   *                        corresponding to the given basis
   */
  static void fromHeterogenerousBasis(const HeterogeneousBasis& basis,
      std::vector<OperationPole*>& operation);

  /**
   * Apply the operator on data.
   *
   * @param[in,out] values    data vector for all full grid points
   *                          (the order is given by IndexVectorRange)
   * @param[in] start         sequence number of the first grid point of the pole
   * @param[in] step          difference of sequence numbers of two subsequent grid points of
   *                          the pole
   * @param[in] count         number of grid points of the pole
   * @param[in] level         level of the full grid
   * @param[in] hasBoundary   whether the full grid has points on the boundary
   */
  void apply(base::DataVector& values, size_t start, size_t step, size_t count,
      level_t level, bool hasBoundary = true) override;

 protected:
  /**
   * Class for the system of linear equations for hierarchising.
   */
  class HierarchisationGeneralSLE : public base::SLE {
   public:
    /**
     * Constructor.
     *
     * @param basis                 1D basis
     * @param dim                   dimensionality
     * @param level                 level of the full grid
     * @param isBasisHierarchical   whether the basis is hierarchical or nodal
     * @param hasBoundary           whether the full grid has points on the boundary
     */
    HierarchisationGeneralSLE(base::Basis<level_t, index_t>& basis, size_t dim, level_t level,
        bool isBasisHierarchical = true, bool hasBoundary = true);

    /**
     * Virtual destructor.
     */
    ~HierarchisationGeneralSLE() override;

    /**
     * Compute an entry of the matrix corresponding to the system of linear equations.
     *
     * @param i   row index
     * @param j   column index
     * @return entry of the matrix
     */
    double getMatrixEntry(size_t i, size_t j) override;

    /**
     * Determine whether an entry of the matrix corresponding to the system of linear equations
     * does not vanish.
     *
     * @param i   row index
     * @param j   column index
     * @return whether the entry of the matrix is non-zero
     */
    bool isMatrixEntryNonZero(size_t i, size_t j) override;

    /**
     * @return whether the basis is hierarchical or nodal
     */
    bool isBasisHierarchical() const;

    /**
     * @param isBasisHierarchical   whether the basis is hierarchical or nodal
     */
    void setIsBasisHierarchical(bool isBasisHierarchical);

    /**
     * @return dimensionality
     */
    size_t getDimension() const override;

    /**
     * @param dim   dimensionality
     */
    void setDimension(size_t dim);

    /**
     * @return level of the full grid
     */
    level_t getLevel() const;

    /**
     * @param level   level of the full grid
     */
    void setLevel(level_t level);

    /**
     * @return whether the full grid has points on the boundary
     */
    bool hasBoundary() const;

    /**
     * @param hasBoundary   whether the full grid has points on the boundary
     */
    void setHasBoundary(bool hasBoundary);

   protected:
    /// 1D basis
    base::Basis<level_t, index_t>& basis;
    /// whether the basis is hierarchical or nodal
    bool isBasisHierarchical_;
    /// dimensionality
    size_t dim;
    /// level of the full grid
    level_t level;
    /// whether the full grid has points on the boundary
    bool hasBoundary_;
  };

  /// system of linear equations for the hierarchising
  HierarchisationGeneralSLE sle;
  /// solver for the system of linear equations
  base::sle_solver::Auto sleSolver;
};

}  // namespace combigrid
}  // namespace sgpp
