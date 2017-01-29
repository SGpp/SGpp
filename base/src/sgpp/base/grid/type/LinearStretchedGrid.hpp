// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARSTRETCHEDGRID_HPP
#define LINEARSTRETCHEDGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with linearstretched base functions
 */
class LinearStretchedGrid : public Grid {
 protected:
  /// grid generator
  StandardGridGenerator generator;
  explicit LinearStretchedGrid(std::istream& istr);

 public:
  /**
   * Constructor LinearStretched Grid without boundaries
   *
   * @param dim the dimension of the grid
   */
  explicit LinearStretchedGrid(size_t dim);

  /**
   * Constructor LinearStretched Grid
   *
   * @param BB the BoundingBox of the grid
   */
  explicit LinearStretchedGrid(Stretching& BB);

  /**
   * Destructor
   */
  ~LinearStretchedGrid() override;

  sgpp::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator& getGenerator() override;
  static Grid* unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp

#endif /* LINEARSTRETCHEDGRID_HPP */
