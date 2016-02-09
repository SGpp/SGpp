// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef GRIDSTENCIL_HPP
#define GRIDSTENCIL_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>


namespace SGPP {
namespace base {

/**
 * grid with linear base functions
 */
class GridStencil : public Grid {
 protected:
  explicit GridStencil(std::istream& istr);

 public:
  typedef std::vector<unsigned int> IndexStencil;

  typedef std::vector<float>      WeightStencil;


  /**
   * Constructor Linear Grid without boundaries
   *
   * @param dim the dimension of the grid
   */
  explicit GridStencil(size_t dim);

  /**
   * Constructor Linear Grid
   *
   * @param BB the BoundingBox of the grid
   */
  explicit GridStencil(BoundingBox& BB);

  /**
   * Destructor
   */
  ~GridStencil() override;

  /**
   * Get the surplus stencil, in fact an array of unsigned ints.
   */
  virtual const IndexStencil&
  getSurplusStencil() const;

  /**
   * Get the neighbor stencil, in fact an array of unsigned ints.
   */
  virtual const IndexStencil&
  getNeighborStencil() const;

  /**
   * Get the weight stencil, in fact an array of floats.
   */
  virtual const WeightStencil&
  getWeightStencil() const;


 protected:
  IndexStencil surplusStencil;
  IndexStencil neighborStencil;
  WeightStencil weightStencil;
};

}  // namespace base
}  // namespace SGPP

#endif /* GRIDSTENCIL_HPP */
