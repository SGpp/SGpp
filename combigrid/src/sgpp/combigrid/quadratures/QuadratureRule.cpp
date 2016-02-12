// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/quadratures/QuadratureRule.hpp>
#include <iostream>
#include <vector>

// to add another stretching-quadrature rule pair it should suffice to
// increase the NUM_QUADS to NUM_QUADS+1, append the corresponding tuple
// at the end of the tuples array and modify the calculate coefficients method
// accordingly
#define NUM_QUADS 4
combigrid::QRtouple tuples[NUM_QUADS] = {{combigrid::EQUIDISTANT, 0},
                                         {combigrid::CHEBYSHEV, 1},
                                         {combigrid::LEGENDRE, 2},
                                         {combigrid::BASU, 3}};

template <typename _Tp>
combigrid::QuadratureRule<_Tp>::QuadratureRule(int max_lvl) {
  if (max_lvl >= 2) {  // 5 pts?
    MAX_LEVELS = std::vector<int>(NUM_QUADS, max_lvl);
  } else {
    COMBIGRID_OUT_WRN(
        "\n"
        "If your max level is not from 2 GREATER,\n"
        "your chances to get a meaningful numerical result will grow FAINTER\n"
        "Hence I take matters in my own HANDS,\n"
        "and try to fix your ARROGANCE \n"
        "by simply setting max_levels to DEFAULT (MAX_LEVELS = 9),\n"
        "so that your wrong numerical result does not become my FAULT.\n"
        "\n"
        "\n!Your hard working CPU! \n",
        __FILE__, __LINE__)

    MAX_LEVELS = std::vector<int>(NUM_QUADS, 9);
  }

  coefficients = std::vector<_Tp**>(NUM_QUADS, NULL);
#pragma omp parallel for schedule(dynamic, 1)

  for (int j = 0; j < NUM_QUADS; j++) {
    coefficients[j] = reinterpret_cast<_Tp**>(malloc(sizeof(_Tp*) * MAX_LEVELS[j]));

    for (int d = 1; d <= MAX_LEVELS[j]; d++) {
      /**
       * For each of the quadrature rules (from 0 to NUM_QUADS)
       *
       * - allocate memory and pre-calculate values of the coefficients from
       *levels 1,2... MAX_LEVELS's
       *
       * */
      QuadratureRule<_Tp>::calculateCoefficients(d, coefficients[j] + d - 1, tuples[j]);
    }
  }
}

template <typename _Tp>
combigrid::QuadratureRule<_Tp>::QuadratureRule(std::vector<int> lvls) {
  coefficients = std::vector<_Tp**>(NUM_QUADS, NULL);

#pragma omp parallel for schedule(dynamic, 1)

  for (int j = 0; j < NUM_QUADS; j++) {
    coefficients[j] = reinterpret_cast<_Tp**>(malloc(sizeof(_Tp*) * MAX_LEVELS[j]));

    for (int d = 1; d <= MAX_LEVELS[j]; d++) {
      /**
       * For each of the quadrature rules (from 0 to NUM_QUADS-1)
       *- allocate memory and pre-calculate values of the coefficients from
       *levels 1,2... MAX_LEVELS's
       *
       * */
      QuadratureRule<_Tp>::calculateCoefficients(d, coefficients[j] + d - 1, tuples[j]);
    }
  }
}

template <typename _Tp>
combigrid::QuadratureRule<_Tp>::~QuadratureRule() {
  for (int j = 0; j < NUM_QUADS; j++) {
    for (int d = 0; d < MAX_LEVELS[j]; d++)
      if (coefficients[j][d] != NULL) {
        free(coefficients[j][d]);
        coefficients[j][d] = NULL;
      }

    free(coefficients[j]);
  }
}

template <typename _Tp>
_Tp combigrid::QuadratureRule<_Tp>::integrate(CombiGrid<_Tp>* grids,
                                              _Tp (*f)(std::vector<double>)) {
  _Tp result = 0.0;
  int dim = grids->getDim();

  /**
   * this flag specifies if the QuadratureRule integrate function should
   * interpolate the abscissas it uses for computing the integral or
   * simply do a straightforward evaluation of the data stored in the underlying
   *full grid.
   *
   *  Interpolate will be true if the user has NOT explicitly specified a
   *function pointer "f"
   *  of a function to evaluate
   */

  bool interpolate = true;
  interpolate = (f != NULL);

  int nr_fullgrids = grids->getNrFullGrids();

/**
 *  evaluate the quadrature in parallel!
 *  coarse level parallelization - each thread from the pool evaluates a full
 * grid quadrature sequentially.
 */
#pragma omp parallel
  {
#pragma omp for reduction(+ : result) schedule(dynamic, 1)

    for (int j = 0; j < nr_fullgrids; j++) {
      if (grids->getFullGrid(j)->isActive()) {
        result += (_Tp)grids->getCoef(j) *
                  quadrature_full_grid(dim, f, grids->getFullGrid(j), interpolate);
      }
    }
  }

  return result;
}

template <typename _Tp>
_Tp combigrid::QuadratureRule<_Tp>::quadrature_full_grid(int dim, _Tp (*f)(std::vector<double>),
                                                         FGridContainer<_Tp>* gridContainer,
                                                         bool interpolate) {
  _Tp result = 0.0f;
  FullGrid<_Tp>* grid = gridContainer->fg();  // obtain a pointer to the fullgrid

  /** select which coefficient matrix shall we choose, depending on the
   * stretching
   *  along any particular direction. **/
  std::vector<int> coef_selector(dim, 0);

  std::vector<Domain1D> domains;
  GridDomain* domain = grid->getDomain();
  // check for fully infinite integration domain.
  std::vector<bool> fully_infinite(dim, false);

  if (domain != NULL) {
    for (int d = 0; d < dim; d++) {
      Domain1D dom_1d = domain->get1DDomain(d);

      if (dom_1d.getMinDomain() == n_INF && dom_1d.getMaxDomain() == p_INF)
        fully_infinite[d] = true;

      coef_selector[d] = stretchingToIdx(dom_1d.getStretchingType(), tuples, NUM_QUADS);

      if (dom_1d.getLevel() > MAX_LEVELS[coef_selector[d]]) {
        COMBIGRID_OUT_ERR(
            " QUADRATURE FAILED: Active full-grid has level greater than the "
            "currently "
            "supported MAX level for the current quadrature! Will return 0.0!",
            __FILE__, __LINE__);
        return (_Tp)0.0;
      }

      domains.push_back(dom_1d);
    }
  }

  /**
   *
   *  I.2 to the function value at point j, superimpose the jacobian of the
   *domain at that point...
   */

  unsigned int num_elem = grid->getNrElements();  // get the total number of grid points

  // precomupte the function values depending on whether or
  // not the user wants to integrate an external function
  // by the function pointer or he/she wants to sample the
  // values stored on the grid...
  std::vector<_Tp> f_values(num_elem, 0.0);
  std::vector<int> indices(dim, 0);

  if (f == NULL) {
    // i.e. if we are not explicitly evaluating a
    // function but rather "probing the values on the grid"
    if (interpolate) {
      std::vector<double> coords(dim, 0);

      for (unsigned int j = 0; j < num_elem; j++) {
        grid->getCoords(j, coords);
        grid->getVectorIndex(j, indices);
        _Tp jacobian_j = 1.0;

        /***
         * Here we add the jacobian in the whole computation... as for infinite
         * and semi-infinite intervals we are
         * not dealing with linear transformations, the jacobian HAS to be
         * computed at each grid point.... Luckily,
         * it is already available as data stored in the Domai1D objects.
         */
        for (unsigned int d = 0; d < indices.size(); d++) {
          jacobian_j *= (_Tp)domains[d].axisJacobian()[indices[d]];
        }

        f_values[j] = grid->eval(coords) * jacobian_j;
      }
    } else {
      for (unsigned int j = 0; j < num_elem; j++) {
        grid->getVectorIndex(j, indices);
        _Tp jacobian_j = 1.0;

        /***
         * Here we add the jacobian in the whole computation... as for infinite
         * and semi-infinite intervals we are
         * not dealing with linear transformations, the jacobian HAS to be
         * computed at each grid point.... Luckily,
         * it is already available as data stored in the Domai1D objects.
         */
        for (unsigned int d = 0; d < indices.size(); d++) {
          jacobian_j *= (_Tp)domains[d].axisJacobian()[indices[d]];
        }

        f_values[j] =
            grid->getElementVector()[j] * jacobian_j;  // simply pick out the value on the grid...
      }
    }
  } else {
    // i.e. if the user has provided an external function
    // for evaluation, perform the evaluation for now..
    std::vector<double> coords(dim, 0);

    for (unsigned int j = 0; j < num_elem; j++) {
      grid->getCoords(j, coords);
      grid->getVectorIndex(j, indices);
      _Tp jacobian_j = 1.0;

      /***
       * Here we add the jacobian in the whole computation... as for infinite
       * and semi-infinite intervals we are
       * not dealing with linear transformations, the jacobian HAS to be
       * computed at each grid point.... Luckily,
       * it is already available as data stored in the Domai1D objects.
       */
      for (unsigned int d = 0; d < indices.size(); d++) {
        jacobian_j *= (_Tp)domains[d].axisJacobian()[indices[d]];
      }

      f_values[j] = f(coords) * jacobian_j;  // simply pick out the value on the grid...
    }
  }

  std::vector<int> levels = gridContainer->getFGLevels();

  /***
   * At this point of the evaluation we already have the functional values
   *evaluated at, hopefully,
   * the correct abscissas, as well as multiplied by the derivative of the
   *underlying coordinate transformation...
   *
   *  Now, all we need to complete the integration is simply multiply the
   *functional values by the weights
   *  and sum the results up...
   */

  for (unsigned int j = 0; j < num_elem; j++) {
    /// obtain the vector indices to calculate the coefficient...
    grid->getVectorIndex(j, indices);
    _Tp W_j = 1;

    for (unsigned int d = 0; d < indices.size(); d++) {
      int idx_d = indices[d];
      int level = levels[d];

      if (idxToStretching(coef_selector[d], tuples, NUM_QUADS) == BASU && fully_infinite[d]) {
        int N_0 = powerOfTwo[level - 1];

        if (idx_d <= N_0)
          W_j *= (coefficients[coef_selector[d]])[level - 2][N_0 - idx_d];
        else
          W_j *= (coefficients[coef_selector[d]])[level - 2][idx_d - N_0];
      } else {
        W_j *= (coefficients[coef_selector[d]])[level - 1][idx_d];
      }
    }

    // we need to do explicit casting here ... or do we?!
    result += ((_Tp)W_j) * f_values[j];
  }

  return result;
}

template <typename _Tp>
int combigrid::QuadratureRule<_Tp>::stretchingToIdx(Stretching str, QRtouple* qrtuples, int size) {
  int idx = 0;

  for (int i = 0; i < size; i++)
    if (str == qrtuples[i].key) {
      idx = i;
      break;
    }

  return idx;
}

template <typename _Tp>
combigrid::Stretching combigrid::QuadratureRule<_Tp>::idxToStretching(int idx, QRtouple* qrtuples,
                                                                      int size) {
  Stretching str = EQUIDISTANT;

  for (int i = 0; i < size; i++)
    if (idx == qrtuples[i].value) {
      str = qrtuples[i].key;
      break;
    }

  return str;
}

template <typename _Tp>
void combigrid::QuadratureRule<_Tp>::calculateCoefficients(int in_level, _Tp** out_coefs,
                                                           QRtouple rule) {
  switch (rule.key) {
    case EQUIDISTANT:  // trapezoid
      TrapezoidalRule<_Tp>::calculateCoefficients(in_level, out_coefs);
      break;

    case CHEBYSHEV:  // clenshaw-curtis
      ClenshawCurtisQuadrature<_Tp>::calculateCoefficients(in_level, out_coefs);
      break;

    case LEGENDRE:  // gauss-patterson
      GaussPattersonQuadrature<_Tp>::calculateCoefficients(in_level, out_coefs);
      break;

    case BASU:  // basu
      BasuQuadrature<_Tp>::calculateCoefficients(in_level, out_coefs);
      break;

    case ATAN:

    // TODO(Christoph): that might not yield good results
    case TAN:

    // TODO(Christoph): that might not yield good results
    case UNKNOWN:  // trapezoidal rule
      TrapezoidalRule<_Tp>::calculateCoefficients(in_level, out_coefs);
      break;

    default:  // trapezoid quadrature
      TrapezoidalRule<_Tp>::calculateCoefficients(in_level, out_coefs);
      break;
  }
}

template class combigrid::QuadratureRule<float>;
template class combigrid::QuadratureRule<double>;

// add more declarations at your hearth's will !!!!

/** CODE FACTOY - ReFACTORY
 *
 *if current level > max_levels...
 "\nAs I wanted to be fast and efficient,\n"
 "I have pre-computed coefficients,\n"
 "only up to level %d inclusive.\n"
 "But I sense you are getting a little abusive"
 "you want me to calculate, "
 "to number crunch and integrate, "
 "on a 1D grid with level far too great\n "
 "for me to anticipate.\n"
 "So please, sit down and re-adjust,\n"
 "your code or grid structure, if you must,\n"
 "so that we could succeed, "
 "the curse of dimensionality to defeat !\n",MAX_LEVELS);
 *
 *
 */
