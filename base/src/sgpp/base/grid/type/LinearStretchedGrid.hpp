// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARSTRETCHEDGRID_HPP
#define LINEARSTRETCHEDGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with linearstretched base functions
 */
class LinearStretchedGrid : public Grid {
 protected:
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

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  std::unique_ptr<GridGenerator> createGridGenerator() override;
  static std::unique_ptr<Grid> unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace SGPP

#endif /* LINEARSTRETCHEDGRID_HPP */
