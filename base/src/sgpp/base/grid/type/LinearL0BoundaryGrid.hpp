// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARBOUNDARYGRID_HPP
#define LINEARBOUNDARYGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/L0BoundaryGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with linear base functions with boundaries
 */
class LinearL0BoundaryGrid : public Grid {
 protected:
  /// grid generator
  L0BoundaryGridGenerator generator;
  explicit LinearL0BoundaryGrid(std::istream& istr);

 public:
  /**
   * Constructor for the Linear Boundary Grid
   *
   * @param dim the dimension of the grid
   */
  explicit LinearL0BoundaryGrid(size_t dim);

  /**
   * Destructor
   */
  ~LinearL0BoundaryGrid() override;

  sgpp::base::GridType getType() override;

  SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp

#endif /* LINEARBOUNDARYGRID_HPP */
