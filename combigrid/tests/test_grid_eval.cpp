/*
 * grid_eval_test.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: hinojosa
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "sgpp/combigrid/utils/combigrid_utils.hpp"
#include "sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
#include "sgpp/combigrid/combischeme/CombiS_CT.hpp"
#include "sgpp/combigrid/combischeme/CombiTS_CT.hpp"
#include "sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp"
#include "sgpp/combigrid/utils/CombigridLevelVector.hpp"

using namespace combigrid;

/*
 * The function to be interpolated
 */
double f_3D(std::vector<double> coords) {
  return 1.0 + (0.25 * (coords[0] - 0.7) * (coords[0] - 0.7) + 2.0)
      + (0.25 * (coords[1] - 0.7) * (coords[1] - 0.7) + 2.0)
      + (0.25 * (coords[2] - 0.7) * (coords[2] - 0.7) + 2.0);
}

double f_2D(std::vector<double> coords) {

  //  return 1.0 + (0.25 * (coords[0] - 0.7) * (coords[0] - 0.7) + 2.0)
  //        + (0.25 * (coords[1] - 0.7) * (coords[1] - 0.7) + 2.0);

  return 4.0 * (coords[0] * coords[0]) * (coords[1] - coords[1] * coords[1]);

}

BOOST_AUTO_TEST_SUITE(testGridEval)

BOOST_AUTO_TEST_CASE( gridEval )
{

  int dim = 2;
  int level = 5;
  std::vector<double> coords(dim, 0.0);

  // Initialize the combischeme and the boundary points flag vector!
  std::vector<int> levels(dim, level);
  AbstractCombiScheme<float>* scheme = new CombiS_CT<float>(levels);
  std::vector<bool> hasbdrypts(dim, true);

  // now initialize a concrete combigrid implementation object
  CombiGrid<float>* grid = new SerialCombiGrid<float>(dim, hasbdrypts);

  // attach the combi scheme to the current combigrid
  grid->attachCombiScheme(scheme);

  // initialize the combigrid data ( levels_vector + coefficients vector) according to the current combischeme implementation
  grid->re_init();

  // populate the combigrid with fullgrids
  grid->createFullGrids();

  //stretch the grid...
  std::vector<double> min(dim, 0.0);
  std::vector<double> max(dim, 1.0);
  AbstractStretchingMaker* stretching = new CombiEquidistantStretching();
  grid->initializeActiveGridsDomain(min, max, stretching);

  /*get the combined subspaces and coefficients*/
  std::vector < std::vector<int> > levels_vect = grid->getLevelsVector();
  std::vector<float> coeffs = grid->getCoefs();

  unsigned int nr_FG = grid->getNrFullGrids();

  /*fill the array with the fullgrid data from the function*/
  for (unsigned int i = 0; i < nr_FG; i++) {
    FullGrid<float>* fgrid = grid->getFullGrid(i)->fg();
    std::vector<float> data = (fgrid->getElementVector()); // get the points vector
    unsigned int size = data.size();

    for (unsigned int j = 0; j < fgrid->getNrElements(); j++) { // should be the same as data.size()
      fgrid->getCoords(j, coords); // working on unit square ...
      fgrid->getElementVector()[j] = f_2D(coords); // evaluate f on the corresponding point.
    }
  }

  /*evaluate the interpolation by combination at a
   * certain coordinate.
   */

  std::vector<double> eval_coordinates(dim, 0.5);
  double f_val = f_2D(eval_coordinates);
  double n_val = grid->eval(eval_coordinates);
  double err = abs(f_val - n_val);

  BOOST_CHECK_CLOSE( n_val, f_val, 1e-10 );

  scheme = grid->detachCombiScheme();
  delete scheme;

  levels.resize(dim, level); // resize levels as it is deleted at the previous step!
  scheme = new CombiTS_CT<float>(levels);
  grid->attachCombiScheme(scheme);

  // now reinitialize the grid...
  grid->re_init();
  grid->createFullGrids();

  /*get the combined subspaces and coefficients*/
  levels_vect = grid->getLevelsVector();
  coeffs = grid->getCoefs();

  /*fill the array with the fullgrid data from the function*/
  for (unsigned int i = 0; i < grid->getNrFullGrids(); i++) {
    FullGrid<float>* fgrid = grid->getFullGrid(i)->fg();
    std::vector<float> data = (fgrid->getElementVector());

    for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
      fgrid->getCoords(j, coords);
      fgrid->getElementVector()[j] = f_2D(coords);
    }
  }

  /*evaluate the interpolation by combination at a
   * certain coordinate.
   */
  std::cout.precision(8);
  f_val = f_2D(eval_coordinates);
  n_val = grid->eval(eval_coordinates);
  err = abs(f_val - n_val);

  BOOST_CHECK_CLOSE( n_val, f_val, 1e-10 );

  scheme = grid->detachCombiScheme();
  delete scheme;
  // empty the combigrid
  // now test the implementation of the TWO Scale Combination Technique scheme...
  levels_vect = grid->getLevelsVector();
  //  std::vector<int> vec(dim,3);
  //  levels_vect.resize(1,vec);
  scheme = new CombiArbitraryScheme<float>(levels_vect);
  grid->attachCombiScheme(scheme);
  // now reinitialize the grid...

  grid->re_init();
  grid->createFullGrids();
  /*get the combined subspaces and coefficients*/
  levels_vect = grid->getLevelsVector();
  coeffs = grid->getCoefs();

  /*fill the array with the fullgrid data from the function*/
  for (unsigned int i = 0; i < grid->getNrFullGrids(); i++) {
    FullGrid<float>* fgrid = grid->getFullGrid(i)->fg();
    std::vector<float> data = (fgrid->getElementVector());

    for (unsigned int j = 0; j < fgrid->getNrElements(); j++) {
      fgrid->getCoords(j, coords);
      fgrid->getElementVector()[j] = f_2D(coords);
    }
  }

  /*evaluate the interpolation by combination at a
   * certain coordinate.
   */
  f_val = f_2D(eval_coordinates);
  n_val = grid->eval(eval_coordinates);
  err = abs(f_val - n_val);

  BOOST_CHECK_CLOSE( n_val, f_val, 1e-10 );
  grid->deleteFullGrids();
  delete scheme;
  delete grid;
}

BOOST_AUTO_TEST_SUITE_END()

