#include <iostream>
// All SG++ headers
//#include "sgpp_base.hpp"

// Or, better!, include only those that are required
#include "base/datatypes/DataVector.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/GridGenerator.hpp"
#include "base/operation/OperationEval.hpp"
#include "base/operation/BaseOpFactory.hpp"

using namespace std;
using namespace sg;
using namespace sg::base;

// function to interpolate
double f(double x0, double x1) {
  return 16.0 * (x0-1)*x0 * (x1-1)*x1;
}

int main() {
  // create a two-dimensional piecewise bi-linear grid
  int dim = 2;
  Grid* grid = Grid::createLinearGrid(dim);
  GridStorage* gridStorage = grid->getStorage();
  cout << "dimensionality:         " << gridStorage->dim() << endl;

  // create regular grid, level 3
  int level = 3;
  GridGenerator* gridGen = grid->createGridGenerator();
  gridGen->regular(level);
  cout << "number of grid points:  " << gridStorage->size() << endl;

  // create coefficient vector
  DataVector alpha(gridStorage->size());
  alpha.setAll(0.0);
  cout << "length of alpha-vector: " << alpha.getSize() << endl;

  // set function values in alpha
  GridIndex* gp;
  for (int i=0; i < gridStorage->size(); i++) {
    gp = gridStorage->get(i);
    alpha[i] = f(gp->abs(0), gp->abs(1));
  }
  cout << alpha.toString() << endl;

  // hierarchize
  op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);
  cout << alpha.toString() << endl;

  // evaluate
  DataVector p(dim);
  p[0] = 0.52;
  p[1] = 0.73;
  OperationEval* opEval = op_factory::createOperationEval(*grid);
  cout << "u(0.52, 0.73) = " << opEval->eval(alpha, p) << endl;

  delete grid;
}
