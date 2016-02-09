#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/GeneralizedBoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/L0BoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/PeriodicGridGenerator.hpp>
#include <sgpp/base/grid/generation/PrewaveletGridGenerator.hpp>
#include <sgpp/base/grid/generation/SquareRootGridGenerator.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>
#include <sgpp/base/grid/generation/StretchedBoundaryGridGenerator.hpp>

using namespace SGPP::base;

BOOST_AUTO_TEST_CASE(testBoundaryGridGenerator) {
  GridStorage *storage = new GridStorage(2);
  BoundaryGridGenerator* gridgen = new BoundaryGridGenerator(storage);
  
  gridgen->regular(2);
  BOOST_CHECK_EQUAL(storage->size(), 21);
  BOOST_CHECK_EQUAL(gridgen->getNumberOfRefinablePoints(), 12);
  BOOST_CHECK_EQUAL(gridgen->getNumberOfRefinablePointsToMaxLevel(1), 8);
  
  storage->emptyStorage();
  
  gridgen->full(2);
  BOOST_CHECK_EQUAL(storage->size(), 25);
  
  //Not implemented
  //storage->emptyStorage();
  
  //gridgen->cliques(3, 1);
  //BOOST_CHECK_EQUAL(storage->size(), 25);
  
  delete storage;
  delete gridgen;
}
