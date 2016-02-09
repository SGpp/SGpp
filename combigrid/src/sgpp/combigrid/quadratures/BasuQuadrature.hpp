/*
 * BasuQuadrature.hpp
 *
 *  Created on: 26 Sep 2014
 *      Author: kenny
 */

#ifndef BASUQUADRATURE_HPP_
#define BASUQUADRATURE_HPP_
#include <sgpp/combigrid/quadratures/AbstractQuadrature.hpp>
#include <sgpp/combigrid/domain/CombiBasuStretching.hpp>

namespace combigrid {

template <typename _Tp>
class BasuQuadrature : public AbstractQuadratureRule<_Tp> {
 public:
  BasuQuadrature(int max_lvl = 9);

  /*Destructor -> Make sure that all memory allocated for the coefs is deleted*/
  virtual ~BasuQuadrature();

  _Tp integrate(CombiGrid<_Tp>* grids, _Tp (*f)(std::vector<double>) = NULL);

  /**
   *  A public static class method that retrieves the BasuQuadrature's
   *coefficients for level in_level, and stores them
   *  in the out_coefs array
   *
   *  @param in_level - an integer specifying the level for which this
   *quadrature's coefficients shall be calculated.
   *  @param out_coefs - a storage array for the newly computed coefficients.
   *
   */
  static void calculateCoefficients(int in_level, _Tp** out_coefs);

 private:
  int MAX_LEVELS;
  /* Container for the pre-computed coefs for levels 1,2... MAX_LEVELS.
   * coefficients[d] is a pointer to a _Tp array
   * of size 2^(d+1) + 1 containing the basu coefficients for a 1-D grid of size
   * 2^(d+1) + 1*/
  _Tp** coefficients;

  /**
   * Do a full grid integration over a single CGF type of grid.
   *
   */
  _Tp basu_full_grid(int dim, _Tp (*f)(std::vector<double>),
                     FGridContainer<_Tp>* gridContainer, bool interpolate);
};
}

#endif /* BASUQUADRATURE_HPP_ */
