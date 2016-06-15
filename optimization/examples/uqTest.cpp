#include "sgpp/base/tools/GaussLegendreQuadRule1D.hpp"
#include "sgpp/optimization/operation/OptimizationOpFactory.hpp"
#include <sgpp_base.hpp>
#include <sgpp/globaldef.hpp>
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

// double bsplineQuadrature(Grid& grid, int d, DataVector& alpha){
//   int p = 3;
//   GridStorage& storage = grid.getStorage();
//   GridIndex* gp;
//   double res = 0.0;
//   std::cout <<  storage.getSize() << std::endl;
//   for (size_t i = 0; i < storage.getSize(); i++) {
//     double temp_res = 0.0;
//     gp = storage.get(i);
//     int level = gp->getLevel(d);
//     int index = gp->getIndex(d);
//     int erster_abschnitt = std::max(-((index - p)/2), 0);
//     int letzter_abschnitt = std::min(p, static_cast<int>((pow(2, level-1) + 1) - ((index+1)/2 ) + 1));

//     std::cout << i << "::" << erster_abschnitt << "::" << letzter_abschnitt << std::endl;
//     for (int j = erster_abschnitt; j <= letzter_abschnitt; j++){
//       temp_res += bsplineIntegral(p, j);
//     }
//     res += pow(2, -level)*temp_res*alpha[i];
//   }
//   return res;
// }


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

double u(double x1, double x2){
  return 4*x1*(1-x1)*4*x2*(1-x2);
}


int main(int argc, char **argv) {
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

  std::cout << bbase->getIntegral(1,1) << std::endl;
  std::cout << bbase->getIntegral(2,3) << std::endl;
  std::cout << bbase->getIntegral(3,5) << std::endl;
  double res = 0.0;
  int in = 0;
  int le = 0;

  std::unique_ptr<Grid> int_grid = Grid::createBsplineGrid(1, 1);
  // std::unique_ptr<Grid> int_grid = Grid::createLinearGrid(1);
  int_grid->getGenerator().regular(1);

  // a = base_test.getIntegral(2,1);
  // double x_1 = 0.0;
  // for (size_t j = 0; j < quadLevel; j++){
  //   x_1 = (1/4.0) * (coordinates[j] + 3.0);
  //   a += base_test.eval(2,3, x_1)*weights[j];
  // }
  // a /= 4;
  // std::cout << a << std::endl;



  GridStorage& intStorage = int_grid->getStorage();
  sgpp::base::DataVector alpha_int(intStorage.getSize());
  std::cout << "asd" << std::endl;
  std::unique_ptr<OperationQuadrature> opQ(sgpp::op_factory::createOperationQuadrature(*int_grid));
  std::cout << "asd" << std::endl;
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

      for (size_t j = 0; j < quadLevel; j++){
        x = h * (coordinates[j] + static_cast<double>(in));
        int_res_comp += base.eval(le, in, x)*f_1d(d, x)*weights[j];
      }
      int_res_comp *= h;

      GridIndex* gp_int;
      for (size_t j = 0; j < intStorage.getSize(); j++) {
        gp_int = intStorage.get(j);
        x = gp_int->getCoord(0);
        alpha_int[j] = base.eval(le, in, x)*f_1d(d, x);
      }
      hierarch->doHierarchisation(alpha_int);

      // GridIndex* gp_1;
      // for (size_t j = 0; j < intStorage.getSize(); j++) {

      //   gp_1 = intStorage.get(j);
      //   int levl = gp_1->getLevel(0);
      //   int ind = gp_1->getIndex(0);
      //   // int_res += pow(2, -levl)*alpha_int[j];
      //   int_res += bsplineQuad(7, levl, ind)*alpha_int[j];
      //   // std::cout << levl << "||" << ind << "||" << bsplineQuad(3, levl, ind) << "||" << alpha_int[j] << std::endl;
      // }
      int_res = opQ->doQuadrature(alpha_int);

      std::cout << d << "::" << int_res << "::" << int_res_comp << std::endl;

      prod_res *= int_res;

    }
    // std:: cout << "prod_res:" << prod_res << std::endl;
    double v_i = alpha[i];
    res += v_i*prod_res;
    // std:: cout << "res:" << res << std::endl;
  }
  std::cout << res << std::endl;

  // // res = 0.0;
  // std::unique_ptr<OperationQuadrature> opQ(sgpp::op_factory::createOperationQuadrature(*grid));
  // res = opQ->doQuadrature(alpha);
  // std::cout << "exact integral value:  " << res << std::endl;

}
