// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBoundaryBasis.hpp>

namespace sgpp {
namespace base {

/**
 * grid with Clenshaw-Curtis linear base functions with boundaries, pentagon cut
 */
class PolyClenshawCurtisBoundaryGrid : public Grid {
 protected:
  explicit PolyClenshawCurtisBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor Linear Truncated Boundary Clenshaw-Curtis Grid
   *
   * @param dim the dimension of the grid
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  explicit PolyClenshawCurtisBoundaryGrid(size_t dim, size_t degree, level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  ~PolyClenshawCurtisBoundaryGrid() override;

  sgpp::base::GridType getType() override;

  SBasis& getBasis() override;

  size_t getDegree() const;

  GridGenerator& getGenerator() override;

  static Grid* unserialize(std::istream& istr);

  void serialize(std::ostream& ostr, int version = SERIALIZATION_VERSION) override;

 protected:
  /// grid generator
  BoundaryGridGenerator generator;
  /// 1 + how much levels the boundary is coarser than the main axes
  level_t boundaryLevel;
  /// degree of the basis
  size_t degree;
  /// basis
  std::unique_ptr<SPolyClenshawCurtisBoundaryBase> basis_;
};

}  // namespace base
}  // namespace sgpp
