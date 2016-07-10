// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/combigrid/CombiGrid.hpp>
#include <sgpp/combigrid/combigrid/SerialCombiGrid.hpp>
#include <sgpp/combigrid/combischeme/AbstractCombiScheme.hpp>
#include <sgpp/combigrid/combischeme/CombiS_CT.hpp>
#include <sgpp/combigrid/domain/CombiChebyshevStretching.hpp>
#include <sgpp/combigrid/domain/CombiEquidistantStretching.hpp>
#include <sgpp/combigrid/domain/CombiLegendreStretching.hpp>
#include <sgpp/combigrid/fullgrid/CombiFullGrid.hpp>
#include <sgpp/combigrid/fullgrid/GridContainer.hpp>
#include <sgpp/combigrid/quadratures/AbstractQuadrature.hpp>
#include <sgpp/combigrid/quadratures/ClenshawCurtisQuadrature.hpp>
#include <sgpp/combigrid/quadratures/GaussPattersonQuadrature.hpp>
#include <sgpp/combigrid/quadratures/TrapezoidalRule.hpp>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

int fcalls = 0;

double f(int dim, std::vector<double> coordinates) {
  fcalls++;
  double result = 1.0;
  if (coordinates.size() != static_cast<size_t>(dim)) {
    std::cout << " Size of coordinates vector and specified problem dimension "
                 "do not match. Aborting \n";
    exit(EXIT_FAILURE);
  }

  for (int d = 0; d < dim; d++) result *= (1. + 1.0 / dim) * pow(coordinates[d], 1. / dim);

  return result;
}

double integrateByRule(combigrid::AbstractQuadratureRule<double>* rule,
                       combigrid::CombiGrid<double>* grid) {
  // fill grid with function values.
  double result = 0.0;
  std::vector<double> coords(grid->getDim(), 0.0);
  int dim = grid->getDim();

  int nr_fg = grid->getNrFullGrids();
  for (int g = 0; g < nr_fg; g++) {
    combigrid::FGridContainer<double>* wrapper = grid->getFullGrid(g);
    if (wrapper->isActive()) {
      combigrid::FullGrid<double>* fg = wrapper->fg();
      int nr_elems = fg->getNrElements();
      for (int e = 0; e < nr_elems; e++) {
        fg->getStandardCoordinates(e, coords);
        fg->getElementVector()[e] = f(dim, coords);
      }
    }
  }

  result = rule->integrate(grid, NULL);

  return result;
}

int main(int argc, char** argv) {
  // initialize the combigrid and setup the combischeme...

  int dim = 5;
  int max_lvl = 7;

  std::vector<bool> hasbdries(dim, true);

  double* errors = reinterpret_cast<double*>(malloc(max_lvl * 3 * sizeof(double)));
  int* Fcalls = reinterpret_cast<int*>(malloc(max_lvl * 3 * sizeof(int)));

  double a = 0.;
  double b = 1.;
  // create some quadratures...
  combigrid::AbstractQuadratureRule<double>* trapz =
      new combigrid::TrapezoidalRule<double>(max_lvl);
  combigrid::AbstractQuadratureRule<double>* clenshaw =
      new combigrid::ClenshawCurtisQuadrature<double>(max_lvl);
  combigrid::AbstractQuadratureRule<double>* gauss =
      new combigrid::GaussPattersonQuadrature<double>(max_lvl);

  combigrid::AbstractStretchingMaker* equidistant = new combigrid::CombiEquidistantStretching();
  combigrid::AbstractStretchingMaker* chebishev = new combigrid::CombiChebyshevStretching();
  combigrid::AbstractStretchingMaker* legendre = new combigrid::CombiLegendreStretching();

  for (int l = 1; l <= max_lvl; l++) {
    std::cout << "Computing errors for max level = " << l << "\n";

    std::vector<int> levels(dim, l);
    combigrid::CombiGrid<double>* grid = new combigrid::SerialCombiGrid<double>(dim, hasbdries);
    combigrid::AbstractCombiScheme<double>* scheme = new combigrid::CombiS_CT<double>(levels);

    grid->attachCombiScheme(scheme);
    grid->re_init();
    grid->createFullGrids();

    /**
    * Integration via trapezoidal rule
    *
    */

    /*
     * setup the domain
     *
     */
    std::vector<double> min(dim, a);
    std::vector<double> max(dim, b);
    fcalls = 0;
    grid->initializeActiveGridsDomain(min, max, equidistant);
    *(errors + (l - 1) * 3) = fabs(integrateByRule(trapz, grid) - 1.0);
    *(Fcalls + (l - 1) * 3) = fcalls;

    fcalls = 0;
    grid->initializeActiveGridsDomain(min, max, chebishev);
    *(errors + (l - 1) * 3 + 1) = fabs(integrateByRule(clenshaw, grid) - 1.0);
    *(Fcalls + (l - 1) * 3 + 1) = fcalls;

    fcalls = 0;
    grid->initializeActiveGridsDomain(min, max, legendre);
    *(errors + (l - 1) * 3 + 2) = fabs(integrateByRule(gauss, grid) - 1.0);
    *(Fcalls + (l - 1) * 3 + 2) = fcalls;

    delete grid;
    delete scheme;
  }

  std::cout << "Integration completed!\n";

  // print the errors table...
  for (int i = 0; i < max_lvl; i++) {
    for (int j = 0; j < 3; j++) {
      std::cout << *(Fcalls + i * 3 + j) << " | " << *(errors + i * 3 + j) << "\t";
    }
    std::cout << "\n";
  }
}
