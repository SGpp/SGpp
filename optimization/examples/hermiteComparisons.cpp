

// or, better, include only the ones needed
#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <iostream>





using sgpp::optimization::test_problems::UnconstrainedTestProblem;
using sgpp::optimization::test_problems::TestScalarFunction;

/**
 * Before starting with the <tt>main</tt> function,
 * the function \f$f\f$, which we want to interpolate, is defined.
 */
double f(double x0, double x1) { return 16.0 * (x0 - 1) * x0 * (x1 - 1) * x1; }

int main() {
 


size_t dim = 2;

/**
* create test function
*/
std::unique_ptr<TestScalarFunction> function;

 function =(std::unique_ptr<TestScalarFunction>(
      new sgpp::optimization::test_problems::RosenbrockObjective(dim)));




 
  


   /**
   * create grid
   */
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearGrid(dim));

  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  std::cout << "dimensionality:         " << gridStorage.getDimension() << std::endl;

  /**
    generate grid
   */
  size_t level = 3;
  grid->getGenerator().regular(level);
  std::cout << "number of grid points:  " << gridStorage.getSize() << std::endl;

  /**
    initialize alpha
   */
  sgpp::base::DataVector alpha(gridStorage.getSize());
  alpha.setAll(0.0);
  std::cout << "length of alpha vector: " << alpha.getSize() << std::endl;
 
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector coordinates(dim);
    gp.getStandardCoordinates(coordinates);
    std::cout<< coordinates[0]<<"   "<<coordinates[1]<<std::endl;
    alpha[i] = function->evalUndisplaced(coordinates);
  }

  std::cout << "alpha before hierarchization: " << alpha.toString() << std::endl;

  /**
   * calculate hierarchized alpha values
   */
  sgpp::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);
  std::cout << "alpha after hierarchization:  " << alpha.toString() << std::endl;

  /**
   * Finally, a second DataVector is created which is used as a point to
   * evaluate the sparse grid function at. An object is obtained which
   * provides an evaluation operation (of type sgpp::base::OperationEvaluation),
   * and the sparse grid interpolant is evaluated at \f$\vec{p}\f$,
   * which is close to (but not exactly at) a grid point.
   */
  sgpp::base::DataVector p(dim);
  p[0] = 0.45;
  p[1] = 0.45;
  std::unique_ptr<sgpp::base::OperationEval> opEval(sgpp::op_factory::createOperationEval(*grid));
  std::cout << "u(0.52, 0.73) = " << opEval->eval(alpha, p) << std::endl;
  std::cout << "u(0.52, 0.73) = " << function->evalUndisplaced(p) << std::endl;
  std::cout << "u(0.52, 0.73) = " << f(p[0],p[1]) << std::endl;
}

/**
 * The example results in the following output:
 * \verbinclude tutorial.output.txt
 * It can be clearly seen that the surpluses decay with a factor of 1/4:
 * On the first level, we obtain 1, on the second 1/4, and on the third
 * 1/16 as surpluses.
 */
