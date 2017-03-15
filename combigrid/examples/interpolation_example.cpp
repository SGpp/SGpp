// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
  * \page example_interpolation_example_cpp Interpolation Example
 * This tutorial shows with two simple example functions in two respectively three dimensions how
 * interpolation with combigrids can be done. These functions are
 * \f[
 *   f_{2D}\colon [0, 1]^2 \to \mathbb{R},\quad
 *   f(x_0, x_1) := 4 x_0^2  (x_1-x_1^2)
 * \f]
 * and
  * \f[
 *   f_{3D}\colon [0, 1]^3 \to \mathbb{R},\quad
 *   f(x_0, x_1, x_2) := 1+ \sum_{i=0}^2 (\frac{(x_i - 0.7)^2}{4} +2)
 * \f].
 */

#include <sgpp/combigrid/combigrid/SerialCombiGrid.hpp>
#include <sgpp/combigrid/combischeme/CombiS_CT.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

/**
 * The exemplary objective function \f$f_{3D}\f$
*/
double f_3D(const std::vector<double> &coords) {
  return 1.0 + (0.25 * (coords[0] - 0.7) * (coords[0] - 0.7) + 2.0) +
         (0.25 * (coords[1] - 0.7) * (coords[1] - 0.7) + 2.0) +
         (0.25 * (coords[2] - 0.7) * (coords[2] - 0.7) + 2.0);
}

/**
 * The exemplary objective function \f$f_{2D} \f$
 */
double f_2D(const std::vector<double> &coords) {
  return 4.0 * (coords[0] * coords[0]) * (coords[1] - coords[1] * coords[1]);
}

/**
 * The main function starts by setting the function pointer func either to the two or three
 * dimensional function. The default dimension is two.
 */
int main() {
  int dim = 2;
  double (*func)(const std::vector<double> &);

  if (dim == 2)
    func = &f_2D;
  else
    func = &f_3D;

  /**
   * Initialize level to default value of five and create two vectors. One for levels and one for
   * coordinates.
   */
  int level = 5;
  std::vector<int> levels(dim, level);
  std::vector<double> coord(dim, 0.0);

  /**
   * Create Combination Scheme which shall be used for combination.
   */
  combigrid::AbstractCombiScheme<double> *scheme = new combigrid::CombiS_CT<double>(levels);

  /**
   * Set bool vector if boundary points are considered.
   */
  std::vector<bool> boundary_points(dim, true);

  /**
   * Set the domain of the computation.
   */
  std::vector<double> min(dim, 0.0);
  std::vector<double> max(dim, 1.0);

  /**
   *  Set the stretching of the domain --> in the simplest case the grid is not stretched.
   */
  combigrid::AbstractStretchingMaker *stretching = new combigrid::CombiEquidistantStretching();

  /**
   * Now the combination grid can be initialized (boundary fixed).
   */
  combigrid::CombiGrid<double> *grid = new combigrid::SerialCombiGrid<double>(dim, boundary_points);

  /**
   * A scheme governing which grids are combined can be attached.
   */
  grid->attachCombiScheme(scheme);

  /**
   *  Compute the combination coefficients for the current combination scheme.
   */
  grid->re_init();

  /**
   * Assign the stretching.
   */
  grid->initializeActiveGridsDomain(min, max, stretching);

  /**
   * Initialize the grid storage, which will be filled with data in the next step.
   */
  grid->createFullGrids();

  /**
   * Get the number of grids and loop over all of them. In the outer loop get the respective full
   * grid. In the inner loop compute the coordinates of the grid point in the domain and fill the
   * full grid with the data at these coordinates.
   */
  unsigned int n_grids = grid->getNrFullGrids();
  for (unsigned int i = 0; i < n_grids; i++) {
    combigrid::FullGrid<double> *fg = grid->getFullGrid(i)->fg();
    for (int j = 0; j < fg->getNrElements(); j++) {
      fg->getStandardCoordinates(j, coord);
      fg->getElementVector()[j] = func(coord);
    }
  }

  /**
   * Default evaluation coordinate is the midpoint. Evaluate there and print the
   * evaluation coordinates and the function value.
   */
  std::vector<double> eval_coord(dim, 0.5);
  std::cout << grid->eval(eval_coord) << '\t' << func(eval_coord) << std::endl;

  /**
   * Print the combination coefficients and the status of all grids.
   */
  std::vector<std::vector<int>> levels_vect = grid->getLevelsVector();
  std::vector<double> coeffs = grid->getCoefs();

  for (unsigned int i = 0; i < coeffs.size(); i++) {
    for (unsigned int j = 0; j < levels_vect[i].size(); j++) {
      std::cout << levels_vect[i][j] << '\t';
    }

    std::cout << "|\t" << coeffs[i] << "\t|\tstatus active: " << grid->getFullGrid(i)->isActive()
              << "\n";
  }

  delete grid;
  delete scheme;
  delete stretching;
}
