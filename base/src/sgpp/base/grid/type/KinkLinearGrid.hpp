// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * grid with kinked linear base functions
 */
class KinkLinearGrid : public Grid {
 protected:
  /// grid generator
  StandardGridGenerator generator;
  explicit KinkLinearGrid(std::istream& istr);

 public:
  /**
   * Constructor kinked linear grid
   *
   * @param dim the dimension of the grid
   */
  explicit KinkLinearGrid(size_t dim);

  /**
   * Destructor
   */
  ~KinkLinearGrid() override;

  sgpp::base::GridType getType() override;

  SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static Grid* unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace sgpp
