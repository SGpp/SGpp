// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/quadratures/BasuQuadrature.hpp>
#include <vector>

template <typename _Tp>
combigrid::BasuQuadrature<_Tp>::BasuQuadrature(int max_lvl) {
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
    BasuQuadrature<_Tp>::calculateCoefficients(d, coefficients + d - 1);
  }
}

template <typename _Tp>
combigrid::BasuQuadrature<_Tp>::~BasuQuadrature() {
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
_Tp combigrid::BasuQuadrature<_Tp>::integrate(CombiGrid<_Tp>* grids,
                                              _Tp (*f)(std::vector<double>)) {
  int dim = grids->getDim();

  bool badstretching = false;

  if (grids->getStretchingType() != BASU) {
    COMBIGRID_OUT_WRN(
        "\n You try to fool me yet AGAIN,\n"
        "you force me count without a PLAN,\n"
        "So ... PLEASE,\n"
        "for my internal PEACE,\n"
        "think about the grid-stretching you INSIST ON,\n"
        "IS it CHEBYSHEV/BASU or EQUIDISTANT?\n"
        "For BASU stretching, I PREFER\n"
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
        if (static_cast<int>(grids->getFullGrid(j)->getMaxLevel()) > MAX_LEVELS) {
          // problem with parallelization - if this error occurs
          // the end result might not be set
          COMBIGRID_OUT_ERR(
              "BASU QUADRATURE FAILED: Active full-grid has level greater than "
              "the currently "
              "supported MAX level. ABORTING",
              __FILE__, __LINE__);
          error_flag++;
        } else {
          result +=
              (_Tp)grids->getCoef(j) * basu_full_grid(dim, f, grids->getFullGrid(j), badstretching);
        }
      }
    }
  }

  if (error_flag > 0) result = 0.0;

  return result;
}

template <typename _Tp>
_Tp combigrid::BasuQuadrature<_Tp>::basu_full_grid(int dim, _Tp (*f)(std::vector<double>),
                                                   FGridContainer<_Tp>* gridContainer,
                                                   bool badstretching) {
  _Tp result = 0.0f;
  FullGrid<_Tp>* grid = gridContainer->fg();  // obtain a pointer to the fullgrid
  std::vector<_Tp> f_values;
  CombiBasuStretching stretching = CombiBasuStretching();
  AbstractQuadratureRule<_Tp>::getGridValues(grid, badstretching, &stretching, &f_values, f);

  /**check if we have a fully infinite interval to integrate over!!! */
  std::vector<Domain1D> domains;
  GridDomain* domain = grid->getDomain();
  // check for fully infinite integration domain.
  std::vector<bool> fully_infinite(dim, false);

  if (domain != NULL) {
    for (int d = 0; d < dim; d++) {
      Domain1D dom = domain->get1DDomain(d);

      if (dom.getMinDomain() == n_INF && dom.getMaxDomain() == p_INF) fully_infinite[d] = true;
    }
  }

  /*
   * At this point of the evaluation we already have the functional values
   *evaluated at, hopefully,
   * the correct abscissas, as well as multiplied by the derivative of the
   *underlying coordinate transformation...
   *
   *  Now, all we need to complete the integration is simply multiply the
   *functional values by the weights
   *  and sum the results up...
   */
  unsigned int num_elem = grid->getNrElements();  // get the total number of grid points
  std::vector<int> indices(dim, 0);
  std::vector<int> levels = gridContainer->getFGLevels();

  for (unsigned int j = 0; j < num_elem; j++) {
    /// obtain the vector indices to calculate the coefficient...
    grid->getVectorIndex(j, indices);
    _Tp W_j = 1;

    for (unsigned int d = 0; d < indices.size(); d++) {
      int idx_d = indices[d];
      int level = levels[d];

      // if along this dimension we have a fully infinite interval...
      if (fully_infinite[d]) {
        int N_0 = powerOfTwo[level - 1];

        if (idx_d <= N_0)
          W_j *= coefficients[level - 2][N_0 - idx_d];
        else
          W_j *= coefficients[level - 2][idx_d - N_0];

      } else {
        W_j *= coefficients[level - 1][idx_d];
      }
    }

    // we need to do explicit casting here ... or do we?!
    result += ((_Tp)W_j) * f_values[j];
  }

  return result;
}

template <class _Tp>
void combigrid::BasuQuadrature<_Tp>::calculateCoefficients(int in_level, _Tp** out_coefs) {
  int N = powerOfTwo[in_level];
  int n = powerOfTwo[in_level - 1];  // is equal to (N-1)/2;
  *out_coefs = reinterpret_cast<_Tp*>(malloc((N + 1) * sizeof(_Tp)));

  for (int s = 0; s <= N; s++) {
    double factor = 2 * M_PI * s / N;
    double w_s = 1.0 / N;

    for (int p = 1; p < n - 1; p++) w_s += (2.0 / N) * cos(factor * p) / (1 - 4 * p * p);

    w_s += (1.0 / N) * cos(factor * (n - 1)) / (1 - 4 * (n - 1) * (n - 1));

    // we need to add this term as to correct for the fact that we are
    // integrating a function in the general
    // form I = Int_0^inf {g(x)} dx instead in the form  Int_0^inf {f(x)e^-x} dx
    double sec_ = 1 / cos(M_PI * s / (2 * N));
    w_s *= sec_ * sec_;

    (*out_coefs)[s] = (_Tp)w_s;
  }

  (*out_coefs)[0] /= (_Tp)2.0;
  (*out_coefs)[N] /= (_Tp)2.0;
}

template class combigrid::BasuQuadrature<float>;
template class combigrid::BasuQuadrature<double>;
