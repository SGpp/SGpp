#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <math.h>
#include "sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
#include "sgpp/combigrid/combischeme/CombiS_CT.hpp"
#include "sgpp/combigrid/combischeme/CombiTS_CT.hpp"
#include "sgpp/combigrid/fullgrid/CombiFullGrid.hpp"
#include "sgpp/combigrid/quadratures/QuadratureRule.hpp"
#include "sgpp/combigrid/quadratures/ClenshawCurtisQuadrature.hpp"
#include "sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp"

using namespace combigrid;

BOOST_AUTO_TEST_SUITE(testQuadratureRule)

double sin_1D(std::vector<double> coords) {
  return sin(coords[0]);
}

double f_2D(std::vector<double> coords) {

  return sin(coords[0]) * sin(coords[1]);
  //  return sin(coords[0]*M_PI)*sin(coords[1]*M_PI)*sin(coords[2]*M_PI);
}

double f_3D(std::vector<double> coords) {
  double res = exp(-coords[0] * coords[0]) * exp(-coords[1] * coords[1])
      * exp(-coords[2] * coords[2]);

  return res;

}

double f_1D(std::vector<double> coords) {

  double res = exp(-coords[0] * coords[0]);
  return res;

}

BOOST_AUTO_TEST_CASE( testQuadratureRule ){
  int dim = 2;
  int level = 6;
  std::vector<bool> hasbdrypts(dim, true);
  std::vector<int> levels(dim, level);

  CombiGrid<double>* grid = new SerialCombiGrid<double>(dim, hasbdrypts);

  //////// STANDARD COMBI TECHNIQUE
//  std::cout << "########################################\n";
//  std::cout << "#### USING  STD COMBI  SCHEME .......###\n";
//  std::cout << "########################################\n";

  AbstractCombiScheme<double>* scheme = new CombiS_CT<double>(levels);

  // attach the combi scheme to the current combigrid
  grid->attachCombiScheme(scheme);
  // initialize the combigrid data ( levels_vector + coefficients vector) according to the current combischeme implementation
  grid->re_init();
  // populate the combigrid with fullgrids
  grid->createFullGrids();

  // print the bdry points...

  FullGrid<double>* fullg = grid->getFullGrid(0)->fg();

  combigrid::AbstractQuadratureRule<double>* trapz = new QuadratureRule<double>(
      level);

  AbstractQuadratureRule<double>* gauss = new QuadratureRule<double>(level);
  AbstractQuadratureRule<double>* clenshaw = new ClenshawCurtisQuadrature<double>(
      level);

  // initialize the different stretchings...
  AbstractStretchingMaker* chebishev_stretching =
      new CombiChebyshevStretching();
  AbstractStretchingMaker* legendre_stretching =
      new CombiLegendreStretching();
  AbstractStretchingMaker* trapz_stretching =
      new CombiEquidistantStretching();

  // initialize the domain..
  std::vector<double> min(dim, 0.0);
  std::vector<double> max(dim, 1.0);
  //std::cout << "Int(sin(x)sin(y)) inside the [0;1]x[0;1] domain\n";
  double an_res = (cos(0) - cos(1)) * (cos(0) - cos(1));
  grid->initializeActiveGridsDomain(min, max, trapz_stretching);
  double result = trapz->integrate(grid, &f_2D);
  //std::cout << "----> Trapz Result : " << result << "\n";

  BOOST_CHECK_CLOSE(result, an_res, 1e-2);

  grid->initializeActiveGridsDomain(min, max, chebishev_stretching);
  result = clenshaw->integrate(grid, &f_2D);
  //std::cout << "----> Clenshaw-Curtis Result : " << result << "\n";

  BOOST_CHECK_CLOSE(result, an_res, 1e-3);

  grid->initializeActiveGridsDomain(min, max, legendre_stretching);
  result = gauss->integrate(grid, &f_2D);
  //std::cout << "----> GAUSS PATTERSON Result : " << result << "\n";
  //std::cout << "----> Analytical result:" << an_res << "\n";
  scheme = grid->detachCombiScheme();
  delete scheme;

  BOOST_CHECK_CLOSE(result, an_res, 1e-3);

//  std::cout << "########################################\n";
//  std::cout << "#### USING  T_S COMBI  SCHEME .......###\n";
//  std::cout << "########################################\n";

  scheme = new CombiTS_CT<double>(levels);
  grid->attachCombiScheme(scheme);
  grid->re_init();
  grid->createFullGrids();

  //std::cout << "Int(sin(x)sin(y)) inside the [0;1]x[0;1] domain\n";
  grid->initializeActiveGridsDomain(min, max, trapz_stretching);
  result = trapz->integrate(grid, &f_2D);
  //std::cout << "----> TRAPZ Result : " << result << "\n";
  BOOST_CHECK_CLOSE(result, an_res, 1e-2);

  grid->initializeActiveGridsDomain(min, max, chebishev_stretching);
  result = clenshaw->integrate(grid, &f_2D);
  //std::cout << "----> CLENSHAW Result : " << result << "\n";
  BOOST_CHECK_CLOSE(result, an_res, 1e-3);

  grid->initializeActiveGridsDomain(min, max, legendre_stretching);
  result = gauss->integrate(grid, &f_2D);
  //std::cout << "----> GAUSS Result: " << result << "\n";
  //std::cout << "----> Analytical result:" << an_res << "\n";
  BOOST_CHECK_CLOSE(result, an_res, 1e-2);

  scheme = grid->detachCombiScheme();
  delete scheme;

//  std::cout << "########################################\n";
//  std::cout << "#### USING  COMBIARBITRARY SCHEME....###\n";
//  std::cout << "########################################\n";

  std::vector<std::vector<int> > levels_vector = grid->getLevelsVector();
  //create the scheme

  min.clear();
  min.resize(dim, -10.0);
  min[1] = -20;
  max.clear();
  max.resize(dim, -4.0);
  max[1] = 5.0;

  scheme = new CombiArbitraryScheme<double>(levels_vector);
  grid->attachCombiScheme(scheme);
  grid->re_init();
  grid->createFullGrids();
  //std::cout << "Int(sin(x)sin(y)) inside the [-10;-4]x[-20;5] domain\n";
  an_res = (cos(-10) - cos(-4)) * (cos(-20) - cos(5));

  grid->initializeActiveGridsDomain(min, max, trapz_stretching);
  std::vector<double> coords(dim, 0.0);
  unsigned int nr_FG = grid->getNrFullGrids();
  FullGrid<double>* fgrid;

  for (unsigned int i = 0; i < nr_FG; i++) {
    if (grid->getFullGrid(i)->isActive()) {
      fgrid = grid->getFullGrid(i)->fg();
      for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
        fgrid->getCoords(j, coords); // working on unit square ...
        fgrid->getElementVector()[j] = f_2D(coords); // evaluate f on the corresponding point.
      }
    }
  }

  result = trapz->integrate(grid, NULL);
  //std::cout << "---->TRAPZ Result : " << result << "\n";
  BOOST_CHECK_CLOSE( result, an_res, 10 );

  grid->initializeActiveGridsDomain(min, max, chebishev_stretching);
  nr_FG = grid->getNrFullGrids();

  for (unsigned int i = 0; i < nr_FG; i++) {
    if (grid->getFullGrid(i)->isActive()) {
      fgrid = grid->getFullGrid(i)->fg();
      for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
        fgrid->getCoords(j, coords); // working on unit square ...
        fgrid->getElementVector()[j] = f_2D(coords); // evaluate f on the corresponding point.
      }
    }
  }

  result = clenshaw->integrate(grid, NULL);
  //std::cout << "---->CLENSHAW Result : " << result << "\n";
  BOOST_CHECK_CLOSE( result, an_res, 1e-1 );

  grid->initializeActiveGridsDomain(min, max, legendre_stretching);
  nr_FG = grid->getNrFullGrids();

  for (unsigned int i = 0; i < nr_FG; i++) {
    if (grid->getFullGrid(i)->isActive()) {
      fgrid = grid->getFullGrid(i)->fg();
      for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
        fgrid->getCoords(j, coords); // working on unit square ...
        fgrid->getElementVector()[j] = f_2D(coords); // evaluate f on the corresponding point.
      }
    }
  }

  result = gauss->integrate(grid, NULL);
  //std::cout << "---->GAUSS Result : " << result << "\n";
  BOOST_CHECK_CLOSE( result, an_res, 1e-1 );

  //std::cout << "----> Analytical result:" << an_res << "\n";

  delete chebishev_stretching;
  delete legendre_stretching;
  delete trapz_stretching;

  delete gauss;
  delete clenshaw;
  delete trapz;

  delete scheme;
  delete grid;

//  //////// STANDARD COMBI TECHNIQUE
//  std::cout << "#########################################################\n";
//  std::cout << "#### 3DIMENSIONS STD  SCHEME.... INFINITE DOMAIN ###\n";
//  std::cout << "#########################################################\n";

  dim = 3;
  level = 9 ;
  hasbdrypts.resize(dim, true);
  coords.resize(dim, 0.0);
  levels.clear();
  levels.resize(dim, level);

  trapz = new QuadratureRule<double>(level);
  gauss = new QuadratureRule<double>(level);
  clenshaw = new ClenshawCurtisQuadrature<double>(level);
  AbstractQuadratureRule<double>* basu = new QuadratureRule<double>(10);

  grid = new SerialCombiGrid<double>(dim, hasbdrypts);
  scheme = new CombiTS_CT<double>(levels);

  grid->attachCombiScheme(scheme);
  grid->init();
  // setup the infinite domains
  grid->createFullGrids();
  nr_FG = grid->getNrFullGrids();

//  std::cout
//      << "Int(exp(-x**2)exp(-y**2) exp(-z**2)) inside the (-inf;0]x[0;inf)x(-inf;inf) domain\n";
  min.clear();
  max.clear();

  min.push_back(n_INF);
  max.push_back(0);
  min.push_back(0);
  max.push_back(p_INF);
  min.push_back(n_INF);
  max.push_back(0);

  trapz_stretching = new CombiEquidistantStretching();
  grid->initializeActiveGridsDomain(min, max, trapz_stretching);

  for (unsigned int i = 0; i < nr_FG; i++) {
    if (grid->getFullGrid(i)->isActive()) {
      fgrid = grid->getFullGrid(i)->fg();
      for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
        fgrid->getCoords(j, coords); // working on unit square ...
        fgrid->getElementVector()[j] = f_3D(coords); // evaluate f on the corresponding point.
      }
    }
  }

  an_res = sqrt(M_PI) * M_PI / 8; // * 0.13940279264033;
  result = trapz->integrate(grid, NULL);
//  std::cout << "---->TRAPZ Result : " << result << " error:"
//      << abs(result - an_res) << "\n";
  BOOST_CHECK_CLOSE( result, an_res, 1e-3 );

  chebishev_stretching = new CombiChebyshevStretching();
  grid->initializeActiveGridsDomain(min, max, chebishev_stretching);

  for (unsigned int i = 0; i < nr_FG; i++) {
    if (grid->getFullGrid(i)->isActive()) {
      fgrid = grid->getFullGrid(i)->fg();
      for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
        fgrid->getCoords(j, coords); // working on unit square ...
        fgrid->getElementVector()[j] = f_3D(coords); // evaluate f on the corresponding point.
      }
    }
  }

  result = clenshaw->integrate(grid, &f_3D);
  //std::cout << "---->Clenshaw Result : " << result << " error:"
     // << abs(result - an_res) << "\n";
  BOOST_CHECK_CLOSE( result, an_res, 1e-5 );

  legendre_stretching = new CombiLegendreStretching();
  grid->initializeActiveGridsDomain(min, max, legendre_stretching);

  for (unsigned int i = 0; i < nr_FG; i++) {
    if (grid->getFullGrid(i)->isActive()) {
      fgrid = grid->getFullGrid(i)->fg();
      for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
        fgrid->getCoords(j, coords); // working on unit square ...
        fgrid->getElementVector()[j] = f_3D(coords); // evaluate f on the corresponding point.
      }
    }
  }

  result = gauss->integrate(grid, NULL);
  //std::cout << "---->Gauss Result : " << result << " error:"
     // << abs(result - an_res) << "\n";
  BOOST_CHECK_CLOSE( result, an_res, 1e-5 );

  //*******************************///

  AbstractStretchingMaker* basu_stretching = new CombiBasuStretching();
  grid->initializeActiveGridsDomain(min, max, basu_stretching);
  nr_FG = grid->getNrFullGrids();

  for (unsigned int i = 0; i < nr_FG; i++) {
    if (grid->getFullGrid(i)->isActive()) {
      fgrid = grid->getFullGrid(i)->fg();
      for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
        fgrid->getCoords(j, coords); // working on unit square ...
        fgrid->getElementVector()[j] = f_3D(coords); // evaluate f on the corresponding point.
      }
    }
  }

  result = basu->integrate(grid, NULL);
//  std::cout << "---->Basu Result : " << result << " error:"
//      << abs(result - an_res) << "\n";
  BOOST_CHECK_CLOSE( result, an_res, 1e-5 );

  //**************************************
  //std::cout << "---->Analytical Result : " << an_res << "\n";

  delete chebishev_stretching;
  delete legendre_stretching;
  delete basu_stretching;
  delete trapz_stretching;

  delete gauss;
  delete clenshaw;
  delete trapz;
  delete basu;

  delete scheme;
  delete grid;
}

BOOST_AUTO_TEST_SUITE_END()
