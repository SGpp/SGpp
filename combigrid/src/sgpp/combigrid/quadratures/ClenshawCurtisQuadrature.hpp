/*
 * ClenshawCurtisQuadrature.hpp
 *
 *  Created on: 31 Jul 2014
 *      Author: kenny
 */

#ifndef CLENSHAWCURTISQUADRATURE_HPP_
#define CLENSHAWCURTISQUADRATURE_HPP_

#include <sgpp/combigrid/quadratures/AbstractQuadrature.hpp>

namespace combigrid {

template <typename _Tp>
class ClenshawCurtisQuadrature : public AbstractQuadratureRule<_Tp> {
 public:
  /**
   * A constructor for the ClenshawCurtisQuadrature class. When the constructor
   *is called,
   * it pre-computes all coefficients for the 1-dimensional Clenshaw-Curtis
   *quadrature up to a
   * specified maximal resolution
   *
   */
  ClenshawCurtisQuadrature(int max_lvl = 9);

  /*Destructor -> Make sure that all memory allocated for the coefs is deleted*/
  virtual ~ClenshawCurtisQuadrature();

  /**
   * Overwrites the inherited "integrate" function from QuadratureRule. This
   *"ClenshawCurtisQuadrature" class encapsulates the
   * implementation of, well, the Clenshaw-Curtis quadrature rule on a
   *combigrid.
   *
   * @param grids - an implementation of the tamplate Combigrid containing the
   *data we wish to integrate on
   * @param f - an optional function pointer can be handed over if the user
   *wishes to apply the numerical quadrature
   * on a specific function "f" . if no such pointer is supplied, or
   *equivalently a NULL value is given as a 3rd argument
   * of the integrate function, the quadrature rule is applied to whatever
   *information is stored on the points of the
   * underlying full grids...
   */

  _Tp integrate(CombiGrid<_Tp> *grids, _Tp (*f)(std::vector<double>) = NULL);

  /**
   *  A public static class method that retrieves the Clenshaw-Curtis'
   *quadrature coefficients for level in_level, and stores them
   *  in the out_coefs array
   *
   *  @param in_level - an integer specifying the level for which this
   *quadrature's coefficients shall be calculated.
   *  @param out_coefs - a storage array for the newly computed coefficients.
   *
   */
  static void calculateCoefficients(int in_level, _Tp **out_coefs);

 private:
  int MAX_LEVELS;

  /* Container for the pre-computed coefs for levels 1,2... MAX_LEVELS.
   * coefficients[d] is a pointer to a _Tp array
   * of size 2^(d+1) + 1 containing the clenshaw-curtis coefficients for a 1-D
   * grid of size 2^(d+1) + 1*/
  _Tp **coefficients;

  /**
   * Do a full grid integration over a single CGF type of grid.
   *
   */
  _Tp clenshaw_curtis_fullgrid(int dim, _Tp (*f)(std::vector<double>),
                               FGridContainer<_Tp> *gridContainer,
                               bool interpolate);
};
}

#endif /* CLENSHAWCURTISQUADRATURE_HPP_ */
