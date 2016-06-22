#include "sgpp/base/tools/GaussLegendreQuadRule1D.hpp"
#include "sgpp/optimization/operation/OptimizationOpFactory.hpp"
#include <sgpp_base.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_optimization.hpp>

#include <random>
#include <cmath>
#include <algorithm>

using sgpp::optimization::OperationMultipleHierarchisation;
using sgpp::base::DataVector;
using sgpp::base::BsplineGrid;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridIndex;
using sgpp::base::GridStorage;
using sgpp::base::OperationQuadrature;
using sgpp::base::OperationQuadratureMC;


class ExampleFunction : public sgpp::optimization::ScalarFunction {
public:
  /**
   * Constructor.
   */
  ExampleFunction() : sgpp::optimization::ScalarFunction(1) {
  }

  /**
   * Evaluates the test function.
   *
   * @param x     point \f$\vec{x} \in [0, 1]^2\f$
   * @return      \f$f(\vec{x})\f$
   */
  double eval(const sgpp::base::DataVector& x) {
    std::unique_ptr<Grid> grid = Grid::createBsplineBoundaryGrid(1, 3);
    size_t level = 3;
    grid->getGenerator().regular(level);
    GridStorage& gridStorage = grid->getStorage();
    sgpp::base::DataVector alpha(gridStorage.getSize());
    for (size_t i = 0; i < gridStorage.getSize(); i++) {
      GridIndex* gp = gridStorage.get(i);
      double eps = gp->getCoord(0);
      alpha[i] = (x[0] - 0.5) * (x[0] - 0.5) * 2 * eps;
    }
    auto hierarch = sgpp::op_factory::createOperationMultipleHierarchisation(*grid);
    hierarch->doHierarchisation(alpha);
    std::unique_ptr<OperationQuadrature> opQ(sgpp::op_factory::createOperationQuadrature(*grid));
    return opQ->doQuadrature(alpha);
    // return (x[0] - 0.5)*(x[0] - 0.5);
    // return std::sin(8.0 * x[0]) + std::sin(7.0 * x[1]);
  }

  /**
   * @param[out] clone pointer to cloned object
   */
  virtual void clone(
                     std::unique_ptr<sgpp::optimization::ScalarFunction>& clone) const {
    clone = std::unique_ptr<sgpp::optimization::ScalarFunction>(
                                                                new ExampleFunction(*this));
  }
};




double cont_genz(int d, double x){
  double a[2];
  double u[2];
  a[0] = 5;
  a[1] = 5;
  u[0] = 0.5;
  u[1] = 0.5;
  return std::exp( - a[d]*std::abs(x - u[d]));
}

void bsplineQuadTest(int p, int level, int index){
  int erster_abschnitt = std::max(-(index-(p+1)/2), 0);
  int letzter_abschnitt = std::min(p, (1 << level) + (p+1)/2 - index - 1);
  std::cout << erster_abschnitt << "::" << letzter_abschnitt << std::endl;

  erster_abschnitt = std::max(-((index - p)/2), 0);
  letzter_abschnitt = std::min(p, static_cast<int>((pow(2, level-1) + 1) - ((index+1)/2 ) + 1));
  std::cout << erster_abschnitt << "::" << letzter_abschnitt << std::endl;

  for(int a = 0; a <= p; a++){
    if((index - (p+1)/2 + a) >= 0 && (index - (p+1)/2 + a) < pow(2,level)) {
      std::cout << a << ",";
    }
  }
  std::cout << std::endl << "------------" << std::endl;
}


double fac(int n){
  int erg = 1;
  for(int i = 2; i <= n; i++){
    erg *= i;
  }
  return erg;
}

double betaFunkt(int p, int q) {
  return (fac(p-1)*fac(q-1))/fac(p+q-1);
}

double f(int dim, double* x, void* clientData){
  double res = 1.0;
  int alpha_1 = 5;
  int beta_1 = 4;
  res *= (1/betaFunkt(alpha_1,beta_1))*pow(x[0], alpha_1 - 1)*pow(1 - x[0], beta_1 - 1);
  int alpha_2 = 3;
  int beta_2 = 2;
  res *= (1/betaFunkt(alpha_2,beta_2))*pow(x[1], alpha_2 - 1)*pow(1 - x[1], beta_2 - 1);
  return res;
}

double f_1d(int d , double x){
  if(d == 0){
      int alpha_1 = 5;
      int beta_1 = 4;
      return (1/betaFunkt(alpha_1,beta_1))*pow(x, alpha_1 - 1)*pow(1 - x, beta_1 - 1);
    }
    else if(d == 1){
      int alpha_2 = 3;
      int beta_2 = 2;
      return (1/betaFunkt(alpha_2,beta_2))*pow(x, alpha_2 - 1)*pow(1 - x, beta_2 - 1);
    }
    else return 0;
}

double const_one(int d, double x){
  return 1.0;
}

double u(double x1, double x2){
  return 4*x1*(1-x1)*4*x2*(1-x2);
}

double parabel_g(double x , double eps){
  return (x - 0.5)*(x - 0.5)*2*eps;
}

void ew_varianz(){
  // create a two-dimensional piecewise bi-linear grid
  int dim = 2;
  // std::unique_ptr<Grid> grid = Grid::createLinearGrid(dim);
  // std::unique_ptr<Grid> grid = Grid::createBsplineGrid(dim, 3);
  std::unique_ptr<Grid> grid = Grid::createPolyGrid(dim, 2);
  GridStorage& gridStorage = grid->getStorage();

  std::cout << "dimensional:        " << gridStorage.getDimension() << std::endl;

  // create regular grid, level 3
  int level = 3;
  grid->getGenerator().regular(level);
  std::cout << "number of grid points: " << gridStorage.getSize() << std::endl;

  // std::cout << bsplineQuadrature(*grid, 0) << std::endl;
  DataVector alpha(gridStorage.getSize());
  GridIndex* gp;
  double p[2];


  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    gp = gridStorage.get(i);
    p[0] = gp->getCoord(0);
    p[1] = gp->getCoord(1);
    alpha[i] = u(p[0], p[1]);
  }
  sgpp::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);
  // sgpp::op_factory::createOperationMultipleHierarchisation(*grid)->doHierarchisation(alpha);

  sgpp::base::DataVector coordinates;
  sgpp::base::DataVector weights;
  sgpp::base::GaussLegendreQuadRule1D gauss;
  size_t quadLevel = 3;
  gauss.getLevelPointsAndWeights(quadLevel, coordinates, weights);

  sgpp::base::SBasis& base = const_cast<sgpp::base::SBasis&>(grid->getBasis());

  std::unique_ptr<sgpp::base::SBsplineClenshawCurtisBase> bbase;
  bbase.reset(new sgpp::base::SBsplineClenshawCurtisBase(3));

  double res = 0.0;
  int in = 0;
  int le = 0;

  std::unique_ptr<Grid> int_grid = Grid::createBsplineBoundaryGrid(1,3);
  int_grid->getGenerator().regular(7);



  GridStorage& intStorage = int_grid->getStorage();
  sgpp::base::DataVector alpha_int(intStorage.getSize());
  std::unique_ptr<OperationQuadrature> opQ(sgpp::op_factory::createOperationQuadrature(*int_grid));
  auto hierarch = sgpp::op_factory::createOperationMultipleHierarchisation(*int_grid);
  // auto hierarch = sgpp::op_factory::createOperationHierarchisation(*int_grid);

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    gp = gridStorage.get(i);
    p[0] = gp->getCoord(0);
    p[1] = gp->getCoord(1);
    double prod_res = 1.0;

    for (size_t d = 0; d < 2; d++){
      in = gp->getIndex(d);
      le = gp->getLevel(d);
      double h = pow(2, -le);
      //quadrature
      double x = 0.0;
      double int_res = 0.0;
      double int_res_comp = 0.0;

      GridIndex* gp_int;
      for (size_t j = 0; j < intStorage.getSize(); j++) {
        gp_int = intStorage.get(j);
        x = gp_int->getCoord(0);
        alpha_int[j] = base.eval(le, in, x)*base.eval(le, in, x)*f_1d(d, x); // Varianz
        // alpha_int[j] = base.eval(le, in, x)*const_one(d, x); // Erwartungswert
      }
      hierarch->doHierarchisation(alpha_int);
      int_res = opQ->doQuadrature(alpha_int);
      std::cout << d << "::" << int_res << "::" << int_res_comp << std::endl;

      prod_res *= int_res;

    }
    std:: cout << "prod_res:" << prod_res << std::endl;
    
    // double v_i = alpha[i]; // Erwartungswert
    double v_i = alpha[i]*alpha[i]; // Varianz
    res += v_i*prod_res;
    std:: cout << "res:" << res << std::endl;
  }
  res -= 0.711111111*0.711111111;
  std::cout << res << std::endl;
}


void optimize(){
  ExampleFunction f;
  // dimension of domain
  const size_t d = f.getNumberOfParameters();
  // B-spline degree
  const size_t p = 3;
  // maximal number of grid points
  const size_t N = 30;
  // adaptivity of grid generation
  const double gamma = 0.95;
  sgpp::base::ModBsplineGrid grid(d, p);
  sgpp::optimization::IterativeGridGeneratorRitterNovak gridGen(f, grid, N, gamma);
  if (!gridGen.generate()) {
    std::cout << "Grid generation failed, exiting.\n";
  }

  // //////////////////////////////////////////////////////////////////////////
  // HIERARCHIZATION
  // //////////////////////////////////////////////////////////////////////////

  std::cout << "Hierarchizing...\n\n";
  sgpp::base::DataVector functionValues(gridGen.getFunctionValues());
  sgpp::base::DataVector coeffs(functionValues.getSize());
  sgpp::optimization::HierarchisationSLE hierSLE(grid);
  sgpp::optimization::sle_solver::Auto sleSolver;

  // solve linear system
  if (!sleSolver.solve(hierSLE, functionValues, coeffs)) {
    std::cout << "Solving failed, exiting.\n";
  }
  // //////////////////////////////////////////////////////////////////////////
  // OPTIMIZATION OF THE SMOOTH INTERPOLANT
  // //////////////////////////////////////////////////////////////////////////

  std::cout << "Optimizing smooth interpolant...\n\n";
  sgpp::optimization::InterpolantScalarFunction ft(grid, coeffs);
  sgpp::optimization::InterpolantScalarFunctionGradient ftGradient(grid, coeffs);
  sgpp::optimization::optimizer::GradientDescent gradientMethod(ft, ftGradient);
  sgpp::base::DataVector x0(d);
  double fX0;
  double ftX0;

  // determine best grid point as starting point
  {
    sgpp::base::GridStorage& gridStorage = grid.getStorage();

    // index of grid point with minimal function value
    size_t x0Index = std::distance(
                       functionValues.getPointer(),
                       std::min_element(functionValues.getPointer(),
                                        functionValues.getPointer() +
                                        functionValues.getSize()));

    for (size_t t = 0; t < d; t++) {
      x0[t] = gridStorage[x0Index]->getCoord(t);
    }

    fX0 = functionValues[x0Index];
    ftX0 = ft.eval(x0);
  }

  std::cout << "x0 = " << x0.toString() << "\n";
  std::cout << "f(x0) = " << fX0 << ", ft(x0) = " << ftX0 << "\n\n";

  gradientMethod.setStartingPoint(x0);
  gradientMethod.optimize();
  const sgpp::base::DataVector& xOpt = gradientMethod.getOptimalPoint();
  const double ftXOpt = gradientMethod.getOptimalValue();
  const double fXOpt = f.eval(xOpt);

  std::cout << "\nxOpt = " << xOpt.toString() << "\n";
  std::cout << "f(xOpt) = " << fXOpt << ", ft(xOpt) = " << ftXOpt << "\n\n";

}


int main(int argc, char **argv) {
  optimize();
  // ExampleFunction f;
  // DataVector x(1);
  // x[0] = 0.8;
  // std::cout << f.eval(x) << std::endl;
  // ew_varianz();
}
