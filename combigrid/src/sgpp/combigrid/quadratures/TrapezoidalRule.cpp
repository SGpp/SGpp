/*
 * TrapezoidalRule.cpp
 *
 *  Created on: 21 Jul 2014
 *      Author: kenny
 */

#include <sgpp/combigrid/quadratures/TrapezoidalRule.hpp>
#include <iostream>


template<typename _Tp>
combigrid::TrapezoidalRule<_Tp>::TrapezoidalRule(int max_lvl) {

  if (max_lvl >= 2) // 5 pts?
    MAX_LEVELS = max_lvl;
  else {
    COMBIGRID_OUT_WRN(
      "\n"
      " If your max level is not from 2 GREATER,\n"
      "your chances to get a meaningful numerical result will grow FAINTER\n"
      "Hence I take matters in my own HANDS,\n"
      "and try to fix your ARROGANCE\n"
      "by simply setting max_levels to DEFAULT (MAX_LEVELS = 10),\n"
      "so that your wrong numerical result does not become my FAULT.\n"
      "\n"
      "\n!Your hard working CPU! \n", __FILE__, __LINE__)

    MAX_LEVELS = 10;
  }

  coefficients = (_Tp**) malloc(sizeof(_Tp*) * MAX_LEVELS);
  #pragma omp parallel for schedule(dynamic,1)

  for ( int d = 1; d <= MAX_LEVELS; d++) {
    /**
     * allocate memory and pre-calculate values of levels 1,2... MAX_LEVELS's coefs
     * we start from level 1, up to level MAX_LEVELS ...
     * */
    TrapezoidalRule<_Tp>::calculateCoefficients(d, coefficients + d - 1);
  }
}

template<typename _Tp>
combigrid::TrapezoidalRule<_Tp>::~TrapezoidalRule() {
  // deallocate memory for trapezoidal rule quadrature...
  // starting from level 1 2... until max_levels..
  for (int d = 0; d < MAX_LEVELS; d++)
    if (coefficients[d] != NULL) {
      free(coefficients[d]);
      coefficients[d] = NULL;
    }

  free(coefficients);

}

template<typename _Tp>
_Tp combigrid::TrapezoidalRule<_Tp>::integrate(CombiGrid<_Tp>* grids,
    _Tp (*f)(std::vector<double>)) {

  _Tp result = 0.0;
  int dim = grids->getDim();
  /**
     * this flag (badstretching) specifies if the trapz_full_grid function should
     * interpolate the abscissas it uses for computing the integral or
     * simply probe the data stored in the underlying full grid.
     *
     *  Interpolate will be true if the combigrids' stretching type IS NOT EQUIDISTANT, as required by the
     *  quadrature rule.
     */
  bool badstretching = false;

  if (grids->getStretchingType() != EQUIDISTANT) {
    COMBIGRID_OUT_WRN("\n You try to fool me yet AGAIN,\n"
                      "you force me count without a PLAN,\n"
                      "So ... PLEASE,\n"
                      "for my internal PEACE,\n"
                      "think about the grid-stretching you INSIST ON,\n"
                      "IS it CHEBYSHEV/ TANGENTIAL or EQUIDISTANT?\n"
                      "For EQUIDISTANT stretching, I PREFER\n"
                      "so that the numerics I could CONQUER!\n"
                      "\n But, hell, I will do whatever you would LIKE,\n"
                      "I will crunch the numbers with DISLIKE,\n"
                      "but all inaccuracies I CREATE,"
                      "are not my fault - it is thy FATE!\n"
                      "!Your hard working CPU!\n", __FILE__, __LINE__)
    badstretching = true;
  }

  // do some point (i.e. domain transformations here...)
  int nr_fullgrids = grids->getNrFullGrids();
  // evaluate the quadrature in parallel!
  // coarse level parallelization - each thread from the pool evaluates a full grid quadrature sequentially.
  int error_flag = 0;
  #pragma omp parallel shared(error_flag)
  {
    #pragma omp for reduction(+:result) schedule(dynamic,1)

    for (int j = 0; j < nr_fullgrids; j++) {

      if (error_flag > 0)
        continue;

      // as break statements are not allowed with OpenMP for loops
      // if error has occured simply skip through all the iterations until you reach the end...
      if (grids->getFullGrid(j)->isActive()) {
        if (grids->getFullGrid(j)->getMaxLevel() >
            MAX_LEVELS) // problem with parallelization - if this error occurs the
          // end result might not be set
        {
          COMBIGRID_OUT_ERR(
            "TRAPEZOIDAL QUADRATURE FAILED: Active full-grid has level greater than the currently "
            "supported MAX level. ABORTING", __FILE__,
            __LINE__);
          error_flag++;

        } else
          result += (_Tp) grids->getCoef(j)
                    * trapz_full_grid(dim, f, grids->getFullGrid(j),
                                      badstretching);
      }
    }
  }

  if (error_flag > 0)
    result = 0.0;

  return result;

}

template<typename _Tp>
_Tp combigrid::TrapezoidalRule<_Tp>::trapz_full_grid(int dim,
    _Tp (*f)(std::vector<double>), FGridContainer<_Tp>* gridContainer,
    bool badstretching) {

  _Tp result = 0.0f;
  FullGrid<_Tp>* grid = gridContainer->fg(); // obtain a pointer to the fullgrid
  std::vector<_Tp> f_values;
  CombiEquidistantStretching stretching = CombiEquidistantStretching();
  AbstractQuadratureRule<_Tp>::getGridValues(grid, badstretching, &stretching,
      f_values, f);

  /***
   * At this point of the evaluation we already have the functional values evaluated at, hopefully,
   * the correct abscissas, as well as multiplied by the derivative of the underlying coordinate transformation...
   *
   *  Now, all we need to complete the integration is simply multiply the functional values by the weights
   *  and sum the results up...
   */
  unsigned int num_elem =
    grid->getNrElements(); // get the total number of grid points
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

    // we need to do explicit casting here ... or do we?!
    result += ((_Tp) W_j) * f_values[j];
  }

  return result;

}

template<typename _Tp>
void combigrid::TrapezoidalRule<_Tp>::calculateCoefficients(
  int in_level, _Tp** out_coefs) {

  int N = powerOfTwo[in_level] + 1;
  *out_coefs = (_Tp*) malloc(N * sizeof(_Tp));

  double base_omega = 2.0;
  base_omega /= (N > 1) ? N - 1.0 : 1.0;

  // we need to halve the first coeff.

  (*out_coefs)[0] = (_Tp) (base_omega / 2.0);

  for (int j = 1; j < N - 1; j++)
    (*out_coefs)[j] = (_Tp) base_omega;

  // we need to halve the last coeff...
  (*out_coefs)[N - 1] = (_Tp) (base_omega / 2.0);

}

template class combigrid::TrapezoidalRule<float>;
template class combigrid::TrapezoidalRule<double>;

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
