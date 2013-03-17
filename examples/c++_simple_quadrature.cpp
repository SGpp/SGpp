#include <iostream>
// All SG++ headers
#include "sgpp_base.hpp"

using namespace std;
using namespace sg;
using namespace sg::base;

// function to interpolate
double f(int dim, double *x, void *clientdata) {
  double res = 1.0;
  for (int i=0; i<dim; i++) {
    res *= 4*x[i]*(1-x[i]);
  }
  return res;
}

int main() {
  // create a two-dimensional piecewise bi-linear grid
  int dim = 2;
  Grid* grid = Grid::createLinearGrid(dim);
  GridStorage* gridStorage = grid->getStorage();
  cout << "dimensionality:        " << gridStorage->dim() << endl;

  // create regular grid, level 3
  int level = 3;
  GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  cout << "number of grid points: " << gridStorage->size() << endl;

  // create coefficient vector
  DataVector alpha(gridStorage->size());
  GridIndex* gp;
  double p[2];
  for (int i=0; i < gridStorage->size(); i++) {
    gp = gridStorage->get(i);
    p[0] = gp->abs(0);
    p[1] = gp->abs(1);
    alpha[i] = f(2, p, NULL);
  }
  op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);

  // direct quadrature
  OperationQuadrature* opQ = op_factory::createOperationQuadrature(*grid);
  double res = opQ->doQuadrature(alpha);
  cout << "exact integral value:  " << res << endl;
  delete opQ;

  // Monte Carlo quadrature using 100000 paths
  OperationQuadratureMC opMC(*grid, 100000);
  res = opMC.doQuadrature(alpha);
  cout << "Monte Carlo value:     " << res << endl;
  res = opMC.doQuadrature(alpha);
  cout << "Monte Carlo value:     " << res << endl;

  // Monte Carlo quadrature of a function
  res = opMC.doQuadratureFunc(f, NULL);
  cout << "MC value:              " << res << endl;

  // Monte Carlo quadrature of error
  res = opMC.doQuadratureL2Error(f, NULL, alpha);
  cout << "MC L2-error:           " << res << endl;

  delete grid;
}
