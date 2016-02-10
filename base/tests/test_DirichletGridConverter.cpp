#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/common/DirichletGridConverter.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

using namespace SGPP::base;


BOOST_AUTO_TEST_CASE(test_DirichletGridConverter){
	/// number of the boundary grid's grid points
	size_t numTotalGridPoints;
	/// number of the inner grid's grid points
	size_t numInnerGridPoints;
	/**
	* array to store the position of i-th inner point in
	* the boundary grid's coefficients
	*/
	size_t* conCoefArray;
	Grid* innerGridExact;

	int level = 3;
	int dimension = 2;
	Grid* linearBoundaryGrid = Grid::createLinearBoundaryGrid(dimension);
	GridStorage* linearBoundaryGridStorageExact = linearBoundaryGrid->getStorage();
	GridGenerator* gridGen = linearBoundaryGrid->createGridGenerator();
    gridGen->regular(level);

	// determine the number of grid points for both grids
	numTotalGridPoints = linearBoundaryGridStorageExact->size();
	numInnerGridPoints = linearBoundaryGridStorageExact->getNumInnerPoints();

	DataVector boundaryGridCoeffs(numTotalGridPoints);
	GridIndex* gp;
	for(size_t i = 0; i < numTotalGridPoints; ++i){
		gp = linearBoundaryGridStorageExact->get(i);
		boundaryGridCoeffs[i] =
				gp->getCoord(0) * (gp->getCoord(0)) + gp->getCoord(1) * (gp->getCoord(1));
	}

	SGPP::op_factory::createOperationHierarchisation(*linearBoundaryGrid)->
			doHierarchisation(boundaryGridCoeffs);

	// allocate the translation array for the coefficients
	conCoefArray = new size_t[numInnerGridPoints];

	// Get the algorithmic dimensions
	std::vector<size_t> BSalgoDims = linearBoundaryGrid->getAlgorithmicDimensions();

	// create new inner Grid, with one grid point
	innerGridExact = new LinearGrid(*(linearBoundaryGrid->getBoundingBox()));

	// Set algorithmic dimensions for inner Grid
	innerGridExact->setAlgorithmicDimensions(BSalgoDims);

	// create new DataVector for storing the inner grid's coefficients
	DataVector innerCoeffsExact(numInnerGridPoints);

	// Iterate through all grid points and filter inner points
	size_t numInner = 0;

	for (size_t i = 0; i < numTotalGridPoints; ++i) {
		GridIndex* curPoint = linearBoundaryGridStorageExact->get(i);
		if (curPoint->isInnerPoint() == true) {
		  // handle coefficients
		  conCoefArray[numInner] = i;
		  innerCoeffsExact.set(numInner, boundaryGridCoeffs.get(i));
		  numInner++;
		  // insert point into inner grid
		  innerGridExact->getStorage()->insert(*curPoint);
		}
	}

	//Test BUILD
	DirichletGridConverter dirichGridConverter;
	Grid* innerGridActual;
	DataVector* innerCoeffsActual;
	dirichGridConverter.buildInnerGridWithCoefs(
			*linearBoundaryGrid, boundaryGridCoeffs, &innerGridActual, &innerCoeffsActual);

	BOOST_CHECK_EQUAL(innerCoeffsActual->getSize(), innerCoeffsExact.getSize());

	for(size_t i = 0; i < innerCoeffsActual->getSize(); ++i){
		BOOST_CHECK_EQUAL(innerCoeffsActual->get(i), innerCoeffsExact.get(i));
	}

	GridIndex* curPointExact;
	GridIndex* curPointActual;
	int innerPointIndex = 0;
	GridStorage* linearBoundaryGridStorageActual = innerGridActual->getStorage();
	for(size_t i = 0; i < numTotalGridPoints; ++i){
		curPointExact = linearBoundaryGridStorageExact->get(i);
		if (curPointExact->isInnerPoint() == true) {
			curPointActual = linearBoundaryGridStorageActual->get(innerPointIndex);
			for (size_t curDim = 0; curDim < curPointExact->dim(); ++curDim){
				BOOST_CHECK_EQUAL(curPointActual->getCoord(curDim), curPointExact->getCoord(curDim));
			}
			innerPointIndex++;
		}
	}

	//TEST REBUILD
	dirichGridConverter.rebuildInnerGridWithCoefs(
			*linearBoundaryGrid, boundaryGridCoeffs, &innerGridActual, &innerCoeffsActual);

	BOOST_CHECK_EQUAL(innerCoeffsActual->getSize(), innerCoeffsExact.getSize());

	for(size_t i = 0; i < innerCoeffsActual->getSize(); ++i){
		BOOST_CHECK_EQUAL(innerCoeffsActual->get(i), innerCoeffsExact.get(i));
	}

	innerPointIndex = 0;
	linearBoundaryGridStorageActual = innerGridActual->getStorage();
	for(size_t i = 0; i < numTotalGridPoints; ++i){
		curPointExact = linearBoundaryGridStorageExact->get(i);
		if (curPointExact->isInnerPoint() == true) {
			curPointActual = linearBoundaryGridStorageActual->get(innerPointIndex);
			for (size_t curDim = 0; curDim < curPointExact->dim(); ++curDim){
				BOOST_CHECK_EQUAL(curPointActual->getCoord(curDim), curPointExact->getCoord(curDim));
			}
			innerPointIndex++;
		}
	}

	//TEST calcInnerCoefs
	dirichGridConverter.calcInnerCoefs(boundaryGridCoeffs, *innerCoeffsActual);
	for(size_t i = 0; i < innerCoeffsActual->getSize(); ++i){
		BOOST_CHECK_EQUAL(innerCoeffsActual->get(i), innerCoeffsExact.get(i));
	}

	//TEST updateBoundaryCoefs
	DataVector boundaryGridCoeffsActual(boundaryGridCoeffs);
	dirichGridConverter.updateBoundaryCoefs(boundaryGridCoeffsActual, *innerCoeffsActual);
	for(size_t i = 0; i < innerCoeffsActual->getSize(); ++i){
		BOOST_CHECK_EQUAL(boundaryGridCoeffsActual.get(i), boundaryGridCoeffs.get(i));
	}

	delete conCoefArray;
	delete innerGridExact;
}
