// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PHIDPHIDOWNBBLINEARBOUNDARY_HPP
#define PHIDPHIDOWNBBLINEARBOUNDARY_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

/**
 * Implementation of sweep operator (): 1D Down for
 * Bilinearform \f$\int_{x} \phi(x) \frac{\partial \phi(x)}{x} dx\f$
 * on linear boundary grids
 */
class PhidPhiDownBBLinearBoundary : public PhidPhiDownBBLinear {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit PhidPhiDownBBLinearBoundary(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~PhidPhiDownBBLinearBoundary();

  /**
   * This operations performs the calculation of down in the direction of dimension <i>dim</i>
   *
   * For level zero it's assumed, that both ansatz-functions do exist: 0,0 and 0,1
   * If one is missing this code might produce some bad errors (segmentation fault, wrong
   * calculation
   * result)
   * So please assure that both functions do exist!
   *
   * @param source sgpp::base::DataVector that contains the gridpoint's coefficients (values from
   * the vector of the laplace operation)
   * @param result sgpp::base::DataVector that contains the result of the down operation
   * @param index a iterator object of the grid
   * @param dim current fixed dimension of the 'execution direction'
   */
  virtual void operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                          grid_iterator& index, size_t dim);
};

}  // namespace finance
}  // namespace sgpp

#endif /* PHIDPHIDOWNBBLINEARBOUNDARY_HPP */
