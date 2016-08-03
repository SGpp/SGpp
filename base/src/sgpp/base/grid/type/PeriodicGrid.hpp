// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PERIODICGRID_HPP
#define PERIODICGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/PeriodicGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with modified linear base functions
 */
class PeriodicGrid : public Grid {
 protected:
  /// grid generator
  PeriodicGridGenerator generator;
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

  sgpp::base::GridType getType() override;

  GridGenerator& getGenerator() override;

  SBasis& getBasis() override;

  static Grid* unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp

#endif /* PERIODICGRID_HPP */
