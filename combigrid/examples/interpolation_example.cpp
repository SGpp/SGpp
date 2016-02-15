// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/combischeme/CombiS_CT.hpp>
#include <sgpp/combigrid/combigrid/SerialCombiGrid.hpp>

#include <iostream>
#include <vector>
#include <algorithm>

/*
 * The 3D function to be interpolated
 */
double f_3D(const std::vector<double> &coords) {
  return 1.0 + (0.25 * (coords[0] - 0.7) * (coords[0] - 0.7) + 2.0) +
         (0.25 * (coords[1] - 0.7) * (coords[1] - 0.7) + 2.0) +
         (0.25 * (coords[2] - 0.7) * (coords[2] - 0.7) + 2.0);
}

/*
 * The 2D function to be interpolated
 */
double f_2D(const std::vector<double> &coords) {
  return 4.0 * (coords[0] * coords[0]) * (coords[1] - coords[1] * coords[1]);
}

int main() {
  // Fix dimension
  int dim = 2;
  double (*func)(const std::vector<double> &);

  if (dim == 2)
    func = &f_2D;
  else
    func = &f_3D;

  // resolution level vector
  int level = 5;
  std::vector<int> levels(dim, level);
  // coordinates
  std::vector<double> coord(dim, 0.0);
  //

  // Create Combination Scheme which shall be used for combination
  combigrid::AbstractCombiScheme<double> *scheme = new combigrid::CombiS_CT<double>(levels);

  // set bool vector if boundary points are considered
  std::vector<bool> boundary_points(dim, true);

  // set the domain of the computation
  std::vector<double> min(dim, 0.0);
  std::vector<double> max(dim, 1.0);
  // set the stretching of the domain --> in the simplest case the grid is not
  // stretched
  combigrid::AbstractStretchingMaker *stretching = new combigrid::CombiEquidistantStretching();

  // now the combination grid can be initialized (boundary fixed)
  combigrid::CombiGrid<double> *grid = new combigrid::SerialCombiGrid<double>(dim, boundary_points);

  // a scheme can be attached, governing which grids are combined
  grid->attachCombiScheme(scheme);

  // compute the combination coefficients for the current combination scheme
  grid->re_init();

  // assign the stretching
  grid->initializeActiveGridsDomain(min, max, stretching);

  // initialize the grid storage, which is filled with data later
  grid->createFullGrids();

  // fill the grids
  // get number of grids
  unsigned int n_grids = grid->getNrFullGrids();

  // loop over all grids
  for (unsigned int i = 0; i < n_grids; i++) {
    // get the respective full grid
    combigrid::FullGrid<double> *fg = grid->getFullGrid(i)->fg();

    for (int j = 0; j < fg->getNrElements(); j++) {
      // compute the position of the grid point in the domain
      fg->getCoords(j, coord);
      // fill the fullgrid with the data at the position
      fg->getElementVector()[j] = func(coord);
    }
  }

  std::vector<double> eval_coord(dim, 0.5);

  std::cout << grid->eval(eval_coord) << '\t' << func(eval_coord) << std::endl;

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
