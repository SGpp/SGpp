/*
 * quadrature_convergence_analysis.cpp
 *
 *  Created on: 12 Jan 2015
 *      Author: kenny
 */

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

using namespace combigrid;

int fcalls = 0;

double f(int dim, std::vector<double> coordinates) {
  fcalls++;
  double result = 1.0;
  if (coordinates.size() != static_cast<size_t>(dim)) {
    std::cout << " Size of coordinates vector and specified problem dimension "
                 "do not match. Aborting \n";
    exit(EXIT_FAILURE);
  }

  for (int d = 0; d < dim; d++)
    result *= (1. + 1.0 / dim) * pow(coordinates[d], 1. / dim);

  return result;
}

double integrateByRule(AbstractQuadratureRule<double>* rule,
                       CombiGrid<double>* grid) {
  // fill grid with function values.
  double result = 0.0;
  std::vector<double> coords(grid->getDim(), 0.0);
  int dim = grid->getDim();

  int nr_fg = grid->getNrFullGrids();
  for (int g = 0; g < nr_fg; g++) {
    FGridContainer<double>* wrapper = grid->getFullGrid(g);
    if (wrapper->isActive()) {
      FullGrid<double>* fg = wrapper->fg();
      int nr_elems = fg->getNrElements();
      for (int e = 0; e < nr_elems; e++) {
        fg->getCoords(e, coords);
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

  double* errors = (double*)malloc(max_lvl * 3 * sizeof(double));
  int* Fcalls = (int*)malloc(max_lvl * 3 * sizeof(int));

  double a = 0.;
  double b = 1.;
  // create some quadratures...
  AbstractQuadratureRule<double>* trapz = new TrapezoidalRule<double>(max_lvl);
  AbstractQuadratureRule<double>* clenshaw =
      new ClenshawCurtisQuadrature<double>(max_lvl);
  AbstractQuadratureRule<double>* gauss =
      new GaussPattersonQuadrature<double>(max_lvl);

  AbstractStretchingMaker* equidistant = new CombiEquidistantStretching();
  AbstractStretchingMaker* chebishev = new CombiChebyshevStretching();
  AbstractStretchingMaker* legendre = new CombiLegendreStretching();

  for (int l = 1; l <= max_lvl; l++) {
    std::cout << "Computing errors for max level = " << l << "\n";

    std::vector<int> levels(dim, l);
    CombiGrid<double>* grid = new SerialCombiGrid<double>(dim, hasbdries);
    AbstractCombiScheme<double>* scheme = new CombiS_CT<double>(levels);

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
      std::cout << *(Fcalls + i * 3 + j) << " | " << *(errors + i * 3 + j)
                << "\t";
    }
    std::cout << "\n";
  }
}
