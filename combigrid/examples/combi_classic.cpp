// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include "combigrid.hpp"

/*
 * The function to be interpolated
 */
double f_3D(std::vector<double> coords) {
  return 1.0 + (0.25 * (coords[0] - 0.7) * (coords[0] - 0.7) + 2.0)
         + (0.25 * (coords[1] - 0.7) * (coords[1] - 0.7) + 2.0)
         + (0.25 * (coords[2] - 0.7) * (coords[2] - 0.7) + 2.0);
}

using namespace combigrid;
int main() {

  int dim = 3;
  int level = 3;
  std::vector<double> coords(dim, 0.0);

  /*initialize combischeme (S_CT, TS_CT, CombiArbitraryScheme)*/
  CombiSchemeBasis* scheme = new S_CT(dim, level);

  /*get the combined subspaces and coefficients*/
  std::vector<std::vector<int> > levels = scheme->getLevels();
  std::vector<double> coeffs = scheme->getCoef();

  /*print the combined subspaces and coefficients*/
  std::cout << "combined subspaces: " << std::endl;

  for (size_t i = 0; i < coeffs.size(); i++) {
    for (size_t j = 0; j < levels[i].size(); j++) {
      std::cout << levels[i][j] << '\t';
    }

    std::cout << "|\t" << coeffs[i] << std::endl;
  }

  /*initialize combigrid (SerialCombiGrid, AdaptiveSerialCombiGrid)*/
  AbstractCombiGrid* grid = new SerialCombiGrid(scheme, true);

  /*create the respective full grids --> allocate the memory*/
  grid->createFullGrids();


  /*fill the array with the fullgrid data from the function*/
  for (int i = 0;
       i < grid->getNrFullGrid();
       i++
      ) {
    FullGridD* fgrid = grid->getFullGrid(i);
    std::vector<double> data = (fgrid->getElementVector());

    for (int j = 0; j < fgrid->getNrElements(); j++) {
      fgrid->getCoords(j, coords);

      fgrid->getElementVector()[j] = f_3D(coords);
    }
  }

  /*evaluate the interpolation by combination at a
   * certain coordinate.
   */
  std::cout.precision(8);
  std::vector<double> eval_coordinates(dim, 0.1);
  std::cout << "function value:\t" << f_3D(eval_coordinates)
            << std::endl;
  std::cout << "combined value:\t" << grid->eval(eval_coordinates)
            << std::endl;

}
