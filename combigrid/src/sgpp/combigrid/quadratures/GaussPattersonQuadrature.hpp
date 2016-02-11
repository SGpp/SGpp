/* ****************************************************************************
* Copyright (C) 2015 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Petar Tzenov

#ifndef GAUSSPATTERSONQUADRATURE_HPP_
#define GAUSSPATTERSONQUADRATURE_HPP_

#include <sgpp/combigrid/quadratures/AbstractQuadrature.hpp>
#include <vector>

namespace combigrid {

template <typename _Tp>
class GaussPattersonQuadrature : public AbstractQuadratureRule<_Tp> {
 public:
  /**
     * A constructor for the current quadrature class. When the constructor is
    *called,
     * it pre-loads all coefficients for the 1-dimensional Gauss-Patterson
    *quadrature up to a
     * specified maximal resolution
     *Due to the high numerical accuracy of the GP quadrature rule, we have
    *precomputed the abscissas and the
     * quadrature coefficients up to a maximum of level 9, i.e. 513 grid points.
     *
     * For an Npoint quadrature formula , where N = 2^level + 1, the polynomial
    *degree of accuracy
     * of the GP Quadrature should be 3*(N/2) + 1.
     *
     **/
  explicit GaussPattersonQuadrature(int max_lvl = 9);

  /*Destructor -> Make sure that all memory allocated for the coefs is deleted*/
  virtual ~GaussPattersonQuadrature();

  /**
   * Overwrites the inherited "integrate" function from QuadratureRule. This
   *"GausPattersonQuadrature" class encapsulates the
   * implementation of, well, the Gaus-Patterson quadrature rule on a combigrid.
   *Polynomial degree of exactness
   * on a 1D grid with n-pts : 3*n+/-1
   *
   * @param grids - an implementation of the tamplate Combigrid containing the
   *data we wish to intergrate on
   * @param f - an optional function pointer can be handed over if the user
   *wishes to apply the numerical quadrature
   * on a specific function "f" . if no such pointer is supplied, or
   *equivalently a NULL value is given as a 3rd argument
   * of the integrate funciton, the quadrature rule is applied to whatever
   *information is stored on the points of the
   * underlying full grids...
   */
  _Tp integrate(CombiGrid<_Tp>* grids, _Tp (*f)(std::vector<double>) = NULL);

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
  static void calculateCoefficients(int in_level, _Tp** out_coefs);

 private:
  int MAX_LEVELS;
  _Tp** coefficients;

  /*dim-dimensional FULL GRID implementation of the GAUSS-PATTERSON quadrature
   *rule.
   *
   */

  _Tp gauss_patterson_fullGrid(int dim, _Tp (*f)(std::vector<double>),
                               FGridContainer<_Tp>* gridContainer,
                               bool interpolate);
};
}  // namespace combigrid

#endif /* GAUSSPATTERSONQUADRATURE_HPP_ */
