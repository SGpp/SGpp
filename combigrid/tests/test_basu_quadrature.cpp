#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "sgpp/combigrid/utils/combigrid_utils.hpp"
#include "sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
#include "sgpp/combigrid/combischeme/CombiS_CT.hpp"
#include "sgpp/combigrid/quadratures/QuadratureRule.hpp"
#include "sgpp/combigrid/domain/CombiBasuStretching.hpp"
#include "sgpp/combigrid/utils/CombigridLevelVector.hpp"

using namespace combigrid;

// test functions
double f_1D_1(std::vector<double> coordinates) {
  return 1.0 / (coordinates[0] + 4.0) * exp(-coordinates[0]);
}

double f_1D_2(std::vector<double> coordinates) {
  double x = coordinates[0];
  return exp(-x * x) * cos(x);
}

double f_1D_3(std::vector<double> coordinates) {

  return exp(-coordinates[0]) * sin(coordinates[0]);

}

double f_1D_4(std::vector<double> coordinates) {

  return exp(coordinates[0]) * sin(coordinates[0]);
}

// analytical solutions of the integral Int(0,inf) {e^-x f(x)}
const double primer_1_solution = 0.206346;
const double primer_2_solution = 1.3803884;

BOOST_AUTO_TEST_SUITE(testBasuQuadrature)

BOOST_AUTO_TEST_CASE( primer123 )
{
  int dim = 1;
  int level = 10;
  std::vector<int> levels(dim, level);
  std::vector<bool> bdries(dim, true);

  /// initialize combigrid
  SerialCombiGrid<double> grid(dim, bdries);
  CombiS_CT<double> scheme(levels);
  grid.attachCombiScheme(&scheme);
  grid.init();
  grid.createFullGrids();
////
  QuadratureRule<double> quad(level);
  CombiBasuStretching stretching;

  /***Primer 1***/
  std::vector<double> min(dim, 0.0);
  std::vector<double> max(dim, p_INF);
  grid.initializeActiveGridsDomain(min, max, &stretching);
  double num_solution_1 = quad.integrate(&grid, &f_1D_1);

  /***Primer 2***/
  double a = n_INF; double b = p_INF;
  min.clear(); min.resize(dim,a); max.clear(); max.resize(dim,b);
  grid.initializeActiveGridsDomain(min,max,&stretching);
  double num_solution_2 = quad.integrate(&grid,&f_1D_2);

  /***Primer 3***/
  a = 0.0; b = p_INF;
  min.clear(); min.resize(dim,a); max.clear(); max.resize(dim,b);
  grid.initializeActiveGridsDomain(min,max,&stretching);
  double primer_3_solution = -0.5*(exp(-b)*(sin(b)+ cos(b)) - exp(-a)*(sin(a)+cos(a)));
  double num_solution_3 = quad.integrate(&grid,&f_1D_3);

  /***Primer 4 ***/
  a = n_INF; b = 0.0;
  min.clear(); min.resize(dim,a); max.clear(); max.resize(dim,b);
  grid.initializeActiveGridsDomain(min,max,&stretching);
  double primer_4_solution = 0.5*(exp(a)*(sin(-a)+ cos(-a)) - exp(b)*(sin(-b)+cos(-b)));
  double num_solution_4 = quad.integrate(&grid,&f_1D_4);

  BOOST_CHECK_CLOSE( primer_1_solution, num_solution_1, 1e-3 );
  BOOST_CHECK_CLOSE( primer_2_solution, num_solution_2, 1e-3 );
  BOOST_CHECK_CLOSE( primer_3_solution, num_solution_3, 1e-3 );
  BOOST_CHECK_CLOSE( primer_4_solution, num_solution_4, 1e-3 );

}



BOOST_AUTO_TEST_SUITE_END()
