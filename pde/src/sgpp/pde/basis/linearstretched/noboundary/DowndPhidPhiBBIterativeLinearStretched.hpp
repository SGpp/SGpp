// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DOWNDPHIDPHIDOWNBBITERATIVELINEARSTRETCHED_HPP
#define DOWNDPHIDPHIDOWNBBITERATIVELINEARSTRETCHED_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * This class is helper class to implement the complete Down
 * of following bilinearform \f$\int_{x} \frac{\partial \phi(x)}{x} \frac{\partial \phi(x)}{x} dx\f$
 * for a given dimension by an iterative algorithms on adaptive
 * Sparse Grids with linear ansatzfunctions without boundaries.
 *
 * This is possible due to the fact that the operator's
 * matrix has only entries on the diagonal.
 *
 * -> the Up/Down can be implemented by iterating over
 * all ansatzfunctions
 */
class DowndPhidPhiBBIterativeLinearStretched {
 private:
  /// Pointer to the grid's storage object
  sgpp::base::GridStorage* storage;

 public:
  /**
   * Constructor
   *
   * @param storage Pointer to the grid's storage object
   */
  explicit DowndPhidPhiBBIterativeLinearStretched(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  ~DowndPhidPhiBBIterativeLinearStretched();

  /**
   * This operations performs the calculation of Down in the direction of dimension <i>dim</i>
   * of following bilinearform: \f$\int_{x} \frac{\partial \phi(x)}{x} \frac{\partial \phi(x)}{x}
   * dx\f$
   *
   * @param alpha sgpp::base::DataVector that contains the gridpoint's coefficients
   * @param result sgpp::base::DataVector that contains the result of the down operation
   * @param dim current fixed dimension of the 'execution direction'
   */
  virtual void operator()(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result,
                          size_t dim);
};
}  // namespace pde
}  // namespace sgpp

#endif /* DOWNDPHIDPHIDOWNBBITERATIVELINEARSTRETCHED_HPP */
