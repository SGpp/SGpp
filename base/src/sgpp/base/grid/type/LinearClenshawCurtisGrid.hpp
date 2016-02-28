// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARCLENSHAWCURTISGRID_HPP
#define LINEARCLENSHAWCURTISGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with Clenshaw-Curtis linear base functions with boundaries, pentagon cut
 */
class LinearClenshawCurtisGrid : public Grid {
 protected:
  explicit LinearClenshawCurtisGrid(std::istream& istr);

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
  explicit LinearClenshawCurtisGrid(size_t dim, level_t boundaryLevel = 1);

  /**
   * Constructor Linear Truncated Boundary Clenshaw-Curtis Grid
   *
   * @param BB the BoundingBox of the grid
   * @param boundaryLevel 1 + how much levels the boundary is coarser than
   *                      the main axes, 0 means one level finer,
   *                      1 means same level,
   *                      2 means one level coarser, etc.
   */
  explicit LinearClenshawCurtisGrid(BoundingBox& BB, level_t boundaryLevel = 1);

  /**
   * Destructor
   */
  ~LinearClenshawCurtisGrid() override;

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);

  void serialize(std::ostream& ostr) override;

 protected:
  /// grid generator
  BoundaryGridGenerator generator;
  /// 1 + how much levels the boundary is coarser than the main axes
  level_t boundaryLevel;
};

}  // namespace base
}  // namespace SGPP

#endif /* LINEARCLENSHAWCURTISGRID_HPP */
