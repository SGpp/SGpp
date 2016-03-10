// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PHIDPHIDOWNBBLINEARSTRETCHED_HPP
#define PHIDPHIDOWNBBLINEARSTRETCHED_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

/**
 * Implementation of sweep operator (): 1D Down for
 * Bilinearform \f$\int_{x} \phi(x) \frac{\partial \phi(x)}{x} dx\f$
 */
class PhidPhiDownBBLinearStretched {
 protected:
  typedef sgpp::base::GridStorage::grid_iterator grid_iterator;

  /// Pointer to the sgpp::base::GridStorage Object
  sgpp::base::GridStorage* storage;
  /// Pointer to the bounding box Obejct
  sgpp::base::Stretching* stretching;

 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit PhidPhiDownBBLinearStretched(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~PhidPhiDownBBLinearStretched();

  /**
   * This operations performs the calculation of down in the direction of dimension <i>dim</i>
   *
   * @param source sgpp::base::DataVector that contains the gridpoint's coefficients (values from
   * the vector of the laplace operation)
   * @param result sgpp::base::DataVector that contains the result of the down operation
   * @param index a iterator object of the grid
   * @param dim current fixed dimension of the 'execution direction'
   */
  virtual void operator()(sgpp::base::DataVector& source, sgpp::base::DataVector& result,
                          grid_iterator& index, size_t dim);

 protected:
  /**
   * recursive function for the calculation of Down with and without Bounding Box support
   * (since the calculations are independent from a bounding box)
   *
   * @param source sgpp::base::DataVector that contains the coefficients of the ansatzfunction
   * @param result sgpp::base::DataVector in which the result of the operation is stored
   * @param index reference to a griditerator object that is used navigate through the grid
   * @param dim the dimension in which the operation is executed
   * @param fl function value on the left boundary
   * @param fr function value on the right boundary
   */
  void rec(sgpp::base::DataVector& source, sgpp::base::DataVector& result, grid_iterator& index,
           size_t dim, double fl, double fr);
};

}  // namespace finance
}  // namespace sgpp

#endif /* PHIDPHIDOWNBBLINEARSTRETCHED_HPP */
