// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TRAPEZOIDAL_QUADRATURE__HPP
#define TRAPEZOIDAL_QUADRATURE__HPP

#include <sgpp/combigrid/quadratures/AbstractQuadrature.hpp>
#include <sgpp/combigrid/fullgrid/CombiFullGrid.hpp>
#include <vector>

namespace combigrid {
template <typename _Tp>
class TrapezoidalRule : public AbstractQuadratureRule<_Tp> {
 public:
  explicit TrapezoidalRule(int max_lvl = 9);

  virtual ~TrapezoidalRule();

  /**
   * Overwrites the inherited "integrate" function from QuadratureRule. This
   *"TrapezoidalRule" class encapsulates the
   * implementation of, well, the Trapezoidal Rule quadrature rule on a
   *combigrid.
   *
   *  Note that calling the constructor initializes the pre-computation of the
   *trapezoidal rule
   *  coeffs on the interval [-1;1] !!! This is done so
   *  to maintain consistency with the quadrature coefficients for the
   *Clenshaw-Curtis and Gauss-Patterson rules.
   *  if you wish to integrate over an arbitrary interval, PLEASE setup a
   *CombiEquidistantStretching class
   *  for the GridDomain, as this class implements all relevant transformations,
   *in effect allowing numerical integration
   *  over arbitrary finite, semi-infinite or infinite intervals...
   *
   *
   * @param grids - an implementation of the Combigrid template containing the
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
  ** pre-computes the trapezoidal rule coeffs on the interval [-1;1]!!! This is
  *done so
  *  to maintain consistency with the quadrature coefficients for the
  *Clenshaw-Curtis and Gauss Patterson rules.
  *  if you wish to integrate over an arbitrary interval, PLEASE setup a
  *CombiEquidistantStretching class
  *  for the GridDomain, as this class implements all relevant transformations,
  *in effect allowing numerical
  *  over arbitrary finite, semi-infinite or infinite intervals....
  *
  *  @param in_level - an integer specifying the level for which this
  *quadrature's coefficients shall be calculated.
  *  @param out_coefs - a storage array for the newly computed coefficients.
  *
  */
  static void calculateCoefficients(int in_level, _Tp** out_coefs);

 private:
  /**
   * Do a full grid trapezoidal rule over a single dim-dimensional FULL GRID.
   *
   */
  _Tp trapz_full_grid(int dim, _Tp (*f)(std::vector<double>), FGridContainer<_Tp>* gridContainer,
                      bool interpolate);

  int MAX_LEVELS;
  /* Container for the pre-computed coefs for levels 1,2... MAX_LEVELS.
   * coefficients[d] is a pointer to a _Tp array
   * of size 2^(d+1) + 1 containing the trapezoidal rule coefficients for a 1-D
   * grid of size 2^(d+1) + 1*/
  _Tp** coefficients;
};
}  // namespace combigrid
#endif
