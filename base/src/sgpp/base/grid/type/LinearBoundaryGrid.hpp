// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARTRUNCATEDBOUNDARYGRID_HPP
#define LINEARTRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with linear base functions with boundaries, pentagon cut
 */
class LinearBoundaryGrid : public Grid {
 protected:
  explicit LinearBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param dim           the dimension of the grid
   * @param boundaryLevel level at which the boundary points should be
   *                      inserted (default = 1: boundary has same level
   *                      as main axes)
   */
  explicit LinearBoundaryGrid(size_t dim, level_t boundaryLevel = 1);

  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param BB the BoundingBox of the grid
   * @param boundaryLevel level at which the boundary points should be
   *                      inserted (default = 1: boundary has same level
   *                      as main axes)
   */
  explicit LinearBoundaryGrid(BoundingBox& BB, level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  ~LinearBoundaryGrid() override;

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator* createGridGenerator() override;

  static Grid* unserialize(std::istream& istr);

  void serialize(std::ostream& ostr) override;

 protected:
  /// level at which the boundary points should be inserted
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace SGPP

#endif /* LINEARTRUNCATEDBOUNDARYGRID_HPP */
