// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARGENERALIZEDTRUNCATEDBOUNDARYGRID_HPP_
#define LINEARGENERALIZEDTRUNCATEDBOUNDARYGRID_HPP_

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/GeneralizedBoundaryGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with linear base functions with boundaries, pentagon cut
 * Generalization of the LinearL0Boundary and LinearBoundary Grids
 * The sparse grid does contain all fullgrids with |l| < a given level, and l_min > l_user
 * For l_user = 0 we obtain the LinearL0BoundaryGrid and for l_user = 1 we obtain the linear truncated boundary grid
 */
class LinearTruncatedBoundaryGrid : public Grid {
 protected:
  /// grid generator
  GeneralizedBoundaryGridGenerator generator;
  explicit LinearTruncatedBoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param dim the dimension of the grid
   */
  explicit LinearTruncatedBoundaryGrid(size_t dim);

  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param BB the BoundingBox of the grid
   */
  explicit LinearTruncatedBoundaryGrid(BoundingBox& BB);

  /**
   * Destructor
   */
  ~LinearTruncatedBoundaryGrid() override;

  sgpp::base::GridType getType() override;

  SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static Grid* unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp

#endif /* LINEARGENERALIZEDTRUNCATEDBOUNDARYGRID_HPP_ */
