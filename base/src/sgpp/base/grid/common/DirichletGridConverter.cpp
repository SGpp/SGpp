// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/DirichletGridConverter.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>

#include <sgpp/base/exception/generation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <cstring>
#include <vector>


namespace sgpp {
namespace base {

DirichletGridConverter::DirichletGridConverter() : numTotalGridPoints(0),
  numInnerGridPoints(0), conCoefArray(nullptr), bFirstTime(true) {
}

DirichletGridConverter::~DirichletGridConverter() {
  if (this->numInnerGridPoints > 0) {
    delete[] this->conCoefArray;
  }
}

void DirichletGridConverter::buildInnerGridWithCoefs(Grid& boundaryGrid,
    DataVector& boundaryCoefs, Grid** innerGrid, DataVector** innerCoefs) {
  if (this->bFirstTime == true) {
    if (boundaryGrid.getType() == base::GridType::LinearL0Boundary
        || boundaryGrid.getType() == base::GridType::LinearBoundary) {
      GridStorage& myGridStorage = boundaryGrid.getStorage();

      // determine the number of grid points for both grids
      this->numTotalGridPoints = myGridStorage.getSize();
      this->numInnerGridPoints = myGridStorage.getNumberOfInnerPoints();

      // std::cout << "Total Points: " << this->numTotalGridPoints << std::endl;
      // std::cout << "Inner Points: " << this->numInnerGridPoints << std::endl;

      // allocate the translation array for the coefficients
      this->conCoefArray = new size_t[this->numInnerGridPoints];

      // Get the algorithmic dimensions
      std::vector<size_t> BSalgoDims = boundaryGrid.getAlgorithmicDimensions();

      // create new inner Grid, with one grid point
      *innerGrid = new LinearGrid(boundaryGrid.getBoundingBox());

      // Set algorithmic dimensions for inner Grid
      (*innerGrid)->setAlgorithmicDimensions(BSalgoDims);

      // create new DataVector for storing the inner grid's coefficients
      *innerCoefs = new DataVector(this->numInnerGridPoints);

      // Iterate through all grid points and filter inner points
      size_t numInner = 0;

      for (size_t i = 0; i < this->numTotalGridPoints; i++) {
        GridPoint& curPoint = myGridStorage[i];

        if (curPoint.isInnerPoint() == true) {
          // handle coefficients
          this->conCoefArray[numInner] = i;
          (*innerCoefs)->set(numInner, boundaryCoefs.get(i));
          numInner++;
          // insert point into inner grid
          (*innerGrid)->getStorage().insert(curPoint);
        }
      }

      // std::string inGrid;
      // *InnerGrid->serialize(inGrid);
      // std::cout << inGrid << std::endl;

      // (*InnerGrid)->getStorage()->recalcLeafProperty();

      this->bFirstTime = false;
    } else if (boundaryGrid.getType() ==
               base::GridType::LinearStretchedBoundary) {
      GridStorage& myGridStorage = boundaryGrid.getStorage();

      // determine the number of grid points for both grids
      this->numTotalGridPoints = myGridStorage.getSize();
      this->numInnerGridPoints = myGridStorage.getNumberOfInnerPoints();

      // std::cout << "Total Points: " << this->numTotalGridPoints << std::endl;
      // std::cout << "Inner Points: " << this->numInnerGridPoints << std::endl;

      // allocate the translation array for the coefficients
      this->conCoefArray = new size_t[this->numInnerGridPoints];

      // Get the algorithmic dimensions
      std::vector<size_t> BSalgoDims = boundaryGrid.getAlgorithmicDimensions();

      // create new inner Grid, with one grid point
      *innerGrid = new LinearStretchedGrid(boundaryGrid.getStretching());

      // Set algorithmic dimensions for inner Grid
      (*innerGrid)->setAlgorithmicDimensions(BSalgoDims);

      // create new DataVector for storing the inner grid's coefficients
      *innerCoefs = new DataVector(this->numInnerGridPoints);

      // Iterate through all grid points and filter inner points
      size_t numInner = 0;

      for (size_t i = 0; i < this->numTotalGridPoints; i++) {
        GridPoint& curPoint = myGridStorage[i];

        if (curPoint.isInnerPoint() == true) {
          // handle coefficients
          this->conCoefArray[numInner] = i;
          (*innerCoefs)->set(numInner, boundaryCoefs.get(i));
          numInner++;
          // insert point into inner grid
          (*innerGrid)->getStorage().insert(curPoint);
        }
      }

      // std::string inGrid;
      // *InnerGrid->serialize(inGrid);
      // std::cout << inGrid << std::endl;

      // (*InnerGrid)->getStorage()->recalcLeafProperty();

      this->bFirstTime = false;
    } else {
      throw generation_exception(
        "DirichletGridConverter : buildInnerGridWithCoefs : "
        "Boundary Grid is from an unsupported grid type!");
    }
  } else {
    throw generation_exception(
      "DirichletGridConverter : buildInnerGridWithCoefs : "
      "This method can only be called once for one instance!");
  }
}

void DirichletGridConverter::rebuildInnerGridWithCoefs(Grid& boundaryGrid,
    DataVector& boundaryCoefs, Grid** innerGrid, DataVector** innerCoefs) {
  if (this->bFirstTime == false) {
    if (boundaryGrid.getType() == base::GridType::LinearL0Boundary
        || boundaryGrid.getType() == base::GridType::LinearBoundary
        || boundaryGrid.getType() == base::GridType::LinearStretchedBoundary) {
      GridStorage& myGridStorage = boundaryGrid.getStorage();

      // determine the number of grid points for both grids
      this->numTotalGridPoints = myGridStorage.getSize();
      this->numInnerGridPoints = myGridStorage.getNumberOfInnerPoints();

      // allocate the translation array for the coefficients
      delete[] this->conCoefArray;
      this->conCoefArray = new size_t[this->numInnerGridPoints];

      // Get the algorithmic dimensions
      std::vector<size_t> BSalgoDims = boundaryGrid.getAlgorithmicDimensions();

      // create new inner Grid, with one grid point
      (*innerGrid)->getStorage().clear();

      // Set algorithmic dimensions for inner Grid
      (*innerGrid)->setAlgorithmicDimensions(BSalgoDims);

      // create new DataVector for storing the inner grid's coefficients
      delete (*innerCoefs);
      *innerCoefs = new DataVector(this->numInnerGridPoints);

      // Iterate through all grid points and filter inner points
      size_t numInner = 0;

      for (size_t i = 0; i < this->numTotalGridPoints; i++) {
        GridPoint& curPoint = myGridStorage[i];

        if (curPoint.isInnerPoint() == true) {
          // handle coefficients
          this->conCoefArray[numInner] = i;
          (*innerCoefs)->set(numInner, boundaryCoefs.get(i));
          numInner++;
          // insert point into inner grid
          (*innerGrid)->getStorage().insert(curPoint);
        }
      }

      // (*InnerGrid)->getStorage()->recalcLeafProperty();
    } else {
      throw generation_exception(
        "DirichletGridConverter : rebuildInnerGridWithCoefs : "
        "Boundary Grid is from an unsupported grid type!");
    }
  } else {
    throw generation_exception(
      "DirichletGridConverter : rebuildInnerGridWithCoefs : "
      "This method can only be called after initial build!");
  }
}

void DirichletGridConverter::calcInnerCoefs(DataVector& BoundaryCoefs,
    DataVector& InnerCoefs) {
  if (this->numInnerGridPoints > 0) {
    for (size_t i = 0; i < this->numInnerGridPoints; i++) {
      InnerCoefs.set(i, BoundaryCoefs.get(this->conCoefArray[i]));
    }
  } else {
    throw generation_exception(
      "DirichletGridConverter : calcInnerCoefs : "
      "The inner grid has no gridpoints! "
      "Check adaptivity calls and settings!");
  }
}

void DirichletGridConverter::updateBoundaryCoefs(DataVector& BoundaryCoefs,
    DataVector& InnerCoefs) {
  if (this->numInnerGridPoints > 0) {
    for (size_t i = 0; i < this->numInnerGridPoints; i++) {
      BoundaryCoefs.set(this->conCoefArray[i], InnerCoefs.get(i));
    }
  } else {
    throw generation_exception(
      "DirichletGridConverter : updateBoundaryCoefs : "
      "The inner grid has no gridpoints! "
      "Check adaptivity calls and settings!");
  }
}

}  // namespace base
}  // namespace sgpp
