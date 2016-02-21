// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARGRIDSTENCIL_HPP
#define LINEARGRIDSTENCIL_HPP

#include <sgpp/base/grid/type/GridStencil.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with linear base functions
 */
class LinearGridStencil : public GridStencil {
 protected:
  /// grid generator
  StandardGridGenerator generator;
  explicit LinearGridStencil(std::istream& istr);

 public:
  /**
   * Constructor Linear Grid without boundaries
   *
   * @param dim the dimension of the grid
   */
  explicit LinearGridStencil(size_t dim);

  /**
   * Constructor Linear Grid
   *
   * @param BB the BoundingBox of the grid
   */
  explicit LinearGridStencil(BoundingBox& BB);

  /**
   * Destructor
   */
  ~LinearGridStencil() override;

  SGPP::base::GridType getType() override;

  const SBasis& getBasis() override;

  GridGenerator& getGenerator() override;

  static std::unique_ptr<Grid> unserialize(std::istream& istr);
};

}  // namespace base
}  // namespace SGPP

#endif /* LINEARGRIDSTENCIL_HPP */
