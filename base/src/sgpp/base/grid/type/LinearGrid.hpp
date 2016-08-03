// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARGRID_HPP
#define LINEARGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with linear base functions
 */
class LinearGrid : public Grid {
 protected:
  /// grid generator
  StandardGridGenerator generator;
  explicit LinearGrid(std::istream& istr);

 public:
  /**
   * Constructor Linear Grid without boundaries
   *
   * @param dim the dimension of the grid
   */
  explicit LinearGrid(size_t dim);

  /**
   * Constructor Linear Grid
   *
   * @param BB the BoundingBox of the grid
   */
  explicit LinearGrid(BoundingBox& BB);

  /**
   * Destructor
   */
  ~LinearGrid() override;

  sgpp::base::GridType getType() override;

  SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static Grid* unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp

#endif /* LINEARGRID_HPP */
