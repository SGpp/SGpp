// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CONVERTLINEARTOPREWAVELET_HPP
#define CONVERTLINEARTOPREWAVELET_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {



/**
 * Class that implements the transformation of a hierarchical linear sparse grid to a
 * hierarchical prewavelet sparse grid. Therefore the ()operator is implemented in order to use
 * the sweep algorithm for the grid traversal. Let the coefficients from the hat basis be
 * \f$ h_{l,i}\f$ and from the prewavelet basis \f$ u_{l,i} \f$. To calculate the surplusses,
 * temp values are needed:
 * \f[
 * (l,i)\neq G_{n}^{1}:t_{l,i}=-\frac{6}{10}u_{l,i\pm1}+t_{l+1,2i}
 * \f]
 * All temp values for levels greater than the maximal level of the grid are set to 0. The actual
 * transformation is defined by the following tridiagonal equation system:
 * \f{eqnarray*}{
 * \frac{16}{10}u_{l,i}+\frac{4}{10}u_{l,i\pm2}&=&h_{l,i}-t_{l+1,2i}+\frac{1}{2}t_{l+1,2(i\pm1)}\\\frac{12}{10}u_{l,1}+\frac{4}{10}u_{l,3}&=&h_{l,1}-t_{l+1,2}+\frac{1}{2}t_{l+1,4}\\\frac{12}{10}u_{l,2^{l}-1}+\frac{4}{10}u_{l,2^{l}-3}&=&h_{l,2^{l}-1}-t_{l+1,2(2^{l}-1)}+\frac{1}{2}t_{l+1,2(2^{l}-2)}
 * \f}
 *
 * For solving these tridiagonal systems, the method described in
 * http://www.nrbook.com/a/bookcpdf/c2-4.pdf was used. Some odd Fileopen® plugin is needed to open that
 * file ... sorry for that! the picture depicts all needed variables in oder to perform the transformation:
 * \image html prewavelets_hierarch.png "This picture shows all involved gridpoints (red crosses) and temp values (green circles) to calculate the new hierarchical coefficients (red arrows) and new temp values (green arrows)."
 */
class ConvertLinearToPrewavelet {
 protected:
  typedef GridStorage::grid_iterator grid_iterator;
  typedef level_t level_type;
  typedef index_t index_type;

  /// the grid object
  GridStorage& storage;
  GridStorage& shadowstorage;

 public:
  /**
   * Constructor, must be bind to a grid
   *
   * An adaptive grid with prewavelet ansatz functions requires for operations
   * using the up-down algorithm shadow points. These shadow points a needed just
   * for data transport, thus they do not have an influence on the final function.
   * Please refer to sgpp::pde::UpDownOneOpDimWithShadow for more information.
     *
   * @param storage the grid storage object of the the grid, on which the hierarchisation should be executed
   * @param shadowstorage shadow points (see detailed description)
   */
  ConvertLinearToPrewavelet(GridStorage& storage, GridStorage& shadowstorage) :
    storage(storage), shadowstorage(shadowstorage) {
  }

  /**
   * Destructor
   */
  ~ConvertLinearToPrewavelet() {
  }

  /**
   * Converts a given linear base to a prewavelet base.
   */
  void operator()(DataVector& source, DataVector& result,
                  grid_iterator& index, size_t dim);
};

}  // namespace base
}  // namespace sgpp

#endif /* CONVERTLINEARTOPREWAVELET_HPP */
