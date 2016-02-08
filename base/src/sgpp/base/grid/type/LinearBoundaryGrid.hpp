// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARTRUNCATEDBOUNDARYGRID_HPP
#define LINEARTRUNCATEDBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with linear base functions with boundaries, pentagon cut
 */
class LinearBoundaryGrid : public Grid {
 protected:
  LinearBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param dim           the dimension of the grid
   * @param boundaryLevel level at which the boundary points should be
   *                      inserted (default = 1: boundary has same level
   *                      as main axes)
   */
  LinearBoundaryGrid(size_t dim, level_t boundaryLevel = 1);

  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param BB the BoundingBox of the grid
   * @param boundaryLevel level at which the boundary points should be
   *                      inserted (default = 1: boundary has same level
   *                      as main axes)
   */
  LinearBoundaryGrid(BoundingBox& BB, level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  virtual ~LinearBoundaryGrid() override;

  virtual SGPP::base::GridType getType() override;

  virtual const SBasis& getBasis() override;

  virtual GridGenerator* createGridGenerator() override;

  static Grid* unserialize(std::istream& istr);

  virtual void serialize(std::ostream& ostr) override;

 protected:
  /// level at which the boundary points should be inserted
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace SGPP

#endif /* LINEARTRUNCATEDBOUNDARYGRID_HPP */
