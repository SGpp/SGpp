// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MODLINEARGRID_HPP
#define MODLINEARGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with modified linear base functions
 */
class ModLinearGrid : public Grid {
 protected:
  /// grid generator
  StandardGridGenerator generator;
  explicit ModLinearGrid(std::istream& istr);

 public:
  /**
   * Constructor modified linear grid
   *
   * @param dim the dimension of the grid
   */
  explicit ModLinearGrid(size_t dim);

  /**
   * Destructor
   */
  ~ModLinearGrid() override;

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace SGPP

#endif /* MODLINEARGRID_HPP */
