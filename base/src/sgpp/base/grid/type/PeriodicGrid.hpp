// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PERIODICGRID_HPP
#define PERIODICGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with modified linear base functions
 */
class PeriodicGrid : public Grid {
 protected:
  explicit PeriodicGrid(std::istream& istr);

 public:
  /**
   * Constructor modified linear grid
   *
   * @param dim the dimension of the grid
   */
  explicit PeriodicGrid(size_t dim);

  /**
   * Destructor
   */
  ~PeriodicGrid() override;

  SGPP::base::GridType getType() override;

  std::unique_ptr<GridGenerator> createGridGenerator() override;

  const SBasis& getBasis() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace SGPP

#endif /* PERIODICGRID_HPP */
