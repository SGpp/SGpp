// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LINEARGRIDSTENCIL_HPP
#define LINEARGRIDSTENCIL_HPP

#include <sgpp/base/grid/type/GridStencil.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with linear base functions
 */
class LinearGridStencil : public GridStencil {
 protected:
  LinearGridStencil(std::istream& istr);

 public:
  /**
   * Constructor Linear Grid without boundaries
   *
   * @param dim the dimension of the grid
   */
  LinearGridStencil(size_t dim);

  /**
   * Constructor Linear Grid
   *
   * @param BB the BoundingBox of the grid
   */
  LinearGridStencil(BoundingBox& BB);

  /**
   * Destructor
   */
  virtual ~LinearGridStencil() override;

  virtual SGPP::base::GridType getType() override;

  virtual const SBasis& getBasis() override;

  virtual GridGenerator* createGridGenerator() override;

  static Grid* unserialize(std::istream& istr);


};

}  // namespace base
}  // namespace SGPP

#endif /* LINEARGRIDSTENCIL_HPP */
