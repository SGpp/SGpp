// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SQUAREROOTGRID_HPP_
#define SQUAREROOTGRID_HPP_
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with linear base functions with boundaries, pentagon cut
 */
class SquareRootGrid : public Grid {
 protected:
  explicit SquareRootGrid(std::istream& istr);

 public:
  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param dim the dimension of the grid
   */
  explicit SquareRootGrid(size_t dim);

  /**
   * Constructor Linear Truncated Boundary Grid
   *
   * @param BB the BoundingBox of the grid
   */
  explicit SquareRootGrid(BoundingBox& BB);

  /**
   * Destructor
   */
  ~SquareRootGrid() override;

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator* createGridGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace SGPP

#endif /* SQUAREROOTGRID_HPP_ */
