/* ****************************************************************************
* Copyright (C) 2015 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Petar Tzenov

#include <sgpp/combigrid/quadratures/ClenshawCurtisQuadrature.hpp>
#include <vector>
#include <sgpp/globaldef.hpp>


template <typename _Tp>
combigrid::ClenshawCurtisQuadrature<_Tp>::ClenshawCurtisQuadrature(
    int max_lvl) {
  if (max_lvl >= 2) {  // 5 pts?
    MAX_LEVELS = max_lvl;
  } else {
    COMBIGRID_OUT_WRN(
        "\n"
        " If your max level is not from 2 GREATER,\n"
        "your chances to get a meaningful numerical result will grow FAINTER\n"
        "Hence I take matters in my own HANDS,\n"
        "and try to fix your ARROGANCE\n"
        "by simply setting max_levels to DEFAULT (MAX_LEVELS = 9),\n"
        "so that your wrong numerical result does not become my FAULT.\n"
        "\n"
        "Best,\nyour hard working CPU! \n",
        __FILE__, __LINE__)

    MAX_LEVELS = 9;
  }

  coefficients = reinterpret_cast<_Tp**>(malloc(sizeof(_Tp*) * MAX_LEVELS));
#pragma omp parallel for schedule(dynamic, 1)

  for (int d = 1; d <= MAX_LEVELS; d++) {
    /**
     * allocate memory and pre-calculate values of levels 1,2... MAX_LEVELS's
     * coefs
     *  we start from level 1, up to level MAX_LEVELS ...
     * */
    ClenshawCurtisQuadrature<_Tp>::calculateCoefficients(d,
                                                         coefficients + d - 1);
  }
}

template <typename _Tp>
combigrid::ClenshawCurtisQuadrature<_Tp>::~ClenshawCurtisQuadrature() {
  // deallocate memory for clenshaw quadrature...
  // starting from level 1 2... until max_levels..
  for (int d = 0; d < MAX_LEVELS; d++)
    if (coefficients[d] != NULL) {
      free(coefficients[d]);
      coefficients[d] = NULL;
    }

  free(coefficients);
}

template <typename _Tp>
_Tp combigrid::ClenshawCurtisQuadrature<_Tp>::integrate(
    CombiGrid<_Tp>* grids, _Tp (*f)(std::vector<double>)) {
  int dim = grids->getDim();

  /**
   * this flag (badstretching) specifies if the clenshaw_fullGridIntegration
   *function should
   * interpolate the abscissas it uses for computing the integral or
   * simply probe the data stored in the underlying full grid.
   *
   *  Interpolate will be true if the combigrids' stretching type IS NOT
   *CHEBYSHEV, as required by the
   *  quadrature rule.
   */
  bool badstretching = false;

  if (grids->getStretchingType() != CHEBYSHEV) {
    COMBIGRID_OUT_WRN(
        "\n You try to fool me yet AGAIN,\n"
        "you force me count without a PLAN,\n"
        "So ... PLEASE,\n"
        "for my internal PEACE,\n"
        "think about the grid-stretching you INSIST ON,\n"
        "IS it CHEBYSHEV/ TANGENTIAL or EQUIDISTANT?\n"
        "For CHEBYSHEV stretching, I PREFER\n"
        "so that the numerics I could CONQUER!\n"
        "\n But, hell, I will do whatever you would LIKE,\n"
        "I will crunch the numbers with DISLIKE,\n"
        "but all inaccuracies I CREATE,"
        "are not my fault - it is thy FATE!\n"
        "!Your hard working CPU!\n",
        __FILE__, __LINE__)
    badstretching = true;
  }

  // do some point (i.e. domain transformations here...)
  int nr_fullgrids = grids->getNrFullGrids();
  // evaluate the quadrature in parallel!
  // coarse level parallelization - each thread from the pool evaluates a full
  // grid quadrature sequentially.

  _Tp result = 0.0;
  int error_flag = 0;

#pragma omp parallel shared(error_flag)
  {
#pragma omp for reduction(+ : result) schedule(dynamic, 1)

    for (int j = 0; j < nr_fullgrids; j++) {
      if (error_flag > 0) continue;

      // as break statements are not allowed with OpenMP for loops
      // if error has occured simply skip through all the iterations until you
      // reach the end...
      if (grids->getFullGrid(j)->isActive()) {
        if (static_cast<int>(grids->getFullGrid(j)->getMaxLevel()) >
            MAX_LEVELS) {
          // problem with parallelization - if this error occurs
          // the
          // end result might not be set
          COMBIGRID_OUT_ERR(
              "CLENSHAW-CURTIS QUADRATURE FAILED: Active full-grid has level "
              "greater than the currently "
              "supported MAX level. ABORTING!",
              __FILE__, __LINE__);
          error_flag++;
        } else {
          result += (_Tp)grids->getCoef(j) *
                    clenshaw_curtis_fullgrid(dim, f, grids->getFullGrid(j),
                                             badstretching);
        }
      }
    }
  }

  if (error_flag > 0) result = 0.0;

  return result;
}

template <typename _Tp>
_Tp combigrid::ClenshawCurtisQuadrature<_Tp>::clenshaw_curtis_fullgrid(
    int dim, _Tp (*f)(std::vector<double>), FGridContainer<_Tp>* gridContainer,
    bool badstretching) {
  _Tp result = 0.0f;
  FullGrid<_Tp>* grid =
      gridContainer->fg();  // obtain a pointer to the fullgrid
  std::vector<_Tp> f_values;
  CombiChebyshevStretching stretching = CombiChebyshevStretching();
  AbstractQuadratureRule<_Tp>::getGridValues(grid, badstretching, &stretching,
                                             &f_values, f);

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
  unsigned int num_elem =
      grid->getNrElements();  // get the total number of grid points
  std::vector<int> indices(dim, 0);
  std::vector<int> levels = gridContainer->getFGLevels();

  for (unsigned int j = 0; j < num_elem; j++) {
    /// obtain the vector indices to calculate the coefficient...
    grid->getVectorIndex(j, indices);
    _Tp W_j = 1;

    for (unsigned int d = 0; d < indices.size(); d++) {
      int idx_d = indices[d];
      int level = levels[d];
      W_j *= coefficients[level - 1][idx_d];
    }

    result += ((_Tp)W_j) * f_values[j];
  }

  return result;
}

template <class _Tp>
void combigrid::ClenshawCurtisQuadrature<_Tp>::calculateCoefficients(
    int in_level, _Tp** out_coefs) {
  int N = powerOfTwo[in_level] + 1;
  int n = powerOfTwo[in_level - 1];  // is equal to (N-1)/2;

  *out_coefs = reinterpret_cast<_Tp*>(malloc(N * sizeof(_Tp)));

  // looking beautiful....
  (*out_coefs[0]) = (_Tp)(1.0 / (N * (N - 2.0)));

  for (int i = 2; i < N; i++) {
    double w_i = 0.0f;
    double factor = 2.0 * M_PI * (i - 1.0) / (N - 1.0);

    for (int j = 1; j < n; j++)
      w_i += 2.0 * cos(factor * j) / (1.0 - 4.0 * j * j);

    w_i += cos(factor * n) / (1.0 - 4.0 * n * n);
    w_i = 2.0 * (1.0 + w_i) / (N - 1.0);
    (*out_coefs)[i - 1] = (_Tp)w_i;
  }

  (*out_coefs)[N - 1] = (_Tp)(1.0 / (N * (N - 2.0)));
}

template class combigrid::ClenshawCurtisQuadrature<float>;
template class combigrid::ClenshawCurtisQuadrature<double>;
