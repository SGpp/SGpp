/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "base/grid/common/DirichletGridConverter.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/type/LinearStretchedGrid.hpp"

#include "base/exception/generation_exception.hpp"

#include <cstring>

namespace sg {
  namespace base {

    DirichletGridConverter::DirichletGridConverter() : numTotalGridPoints(0), numInnerGridPoints(0), conCoefArray(NULL), bFirstTime(true) {
    }

    DirichletGridConverter::~DirichletGridConverter() {
      if (this->numInnerGridPoints > 0) {
        delete[] this->conCoefArray;
      }
    }

    void DirichletGridConverter::buildInnerGridWithCoefs(Grid& BoundaryGrid, DataVector& BoundaryCoefs, Grid** InnerGrid, DataVector** InnerCoefs) {
      if (this->bFirstTime == true) {
        if (strcmp(BoundaryGrid.getType(), "linearBoundary") == 0 || strcmp(BoundaryGrid.getType(), "linearTrapezoidBoundary") == 0) {
          GridStorage* myGridStorage = BoundaryGrid.getStorage();

          // determine the number of grid points for both grids
          this->numTotalGridPoints = myGridStorage->size();
          this->numInnerGridPoints = myGridStorage->getNumInnerPoints();

          //std::cout << "Total Points: " << this->numTotalGridPoints << std::endl;
          //std::cout << "Inner Points: " << this->numInnerGridPoints << std::endl;

          // allocate the translation array for the coefficients
          this->conCoefArray = new size_t[this->numInnerGridPoints];

          // Get the algorithmic dimensions
          std::vector<size_t> BSalgoDims = BoundaryGrid.getAlgorithmicDimensions();

          // create new inner Grid, with one grid point
          *InnerGrid = new LinearGrid(*BoundaryGrid.getBoundingBox());

          // Set algorithmic dimensions for inner Grid
          (*InnerGrid)->setAlgorithmicDimensions(BSalgoDims);

          // create new DataVector for storing the inner grid's coefficients
          *InnerCoefs = new DataVector(this->numInnerGridPoints);

          // Iterate through all grid points and filter inner points
          size_t numInner = 0;

          for (size_t i = 0; i < this->numTotalGridPoints; i++) {
            GridIndex* curPoint = (*myGridStorage)[i];

            if (curPoint->isInnerPoint() == true) {
              // handle coefficients
              this->conCoefArray[numInner] = i;
              (*InnerCoefs)->set(numInner, BoundaryCoefs.get(i));
              numInner++;
              // insert point into inner grid
              (*InnerGrid)->getStorage()->insert(*curPoint);
            }
          }

          //std::string inGrid;
          //*InnerGrid->serialize(inGrid);
          //std::cout << inGrid << std::endl;

          //(*InnerGrid)->getStorage()->recalcLeafProperty();

          this->bFirstTime = false;
        } else if (strcmp(BoundaryGrid.getType(), "linearStretchedTrapezoidBoundary") == 0) {
          GridStorage* myGridStorage = BoundaryGrid.getStorage();

          // determine the number of grid points for both grids
          this->numTotalGridPoints = myGridStorage->size();
          this->numInnerGridPoints = myGridStorage->getNumInnerPoints();

          //std::cout << "Total Points: " << this->numTotalGridPoints << std::endl;
          //std::cout << "Inner Points: " << this->numInnerGridPoints << std::endl;

          // allocate the translation array for the coefficients
          this->conCoefArray = new size_t[this->numInnerGridPoints];

          // Get the algorithmic dimensions
          std::vector<size_t> BSalgoDims = BoundaryGrid.getAlgorithmicDimensions();

          // create new inner Grid, with one grid point
          *InnerGrid = new LinearStretchedGrid(*BoundaryGrid.getStretching());

          // Set algorithmic dimensions for inner Grid
          (*InnerGrid)->setAlgorithmicDimensions(BSalgoDims);

          // create new DataVector for storing the inner grid's coefficients
          *InnerCoefs = new DataVector(this->numInnerGridPoints);

          // Iterate through all grid points and filter inner points
          size_t numInner = 0;

          for (size_t i = 0; i < this->numTotalGridPoints; i++) {
            GridIndex* curPoint = (*myGridStorage)[i];

            if (curPoint->isInnerPoint() == true) {
              // handle coefficients
              this->conCoefArray[numInner] = i;
              (*InnerCoefs)->set(numInner, BoundaryCoefs.get(i));
              numInner++;
              // insert point into inner grid
              (*InnerGrid)->getStorage()->insert(*curPoint);
            }
          }

          //std::string inGrid;
          //*InnerGrid->serialize(inGrid);
          //std::cout << inGrid << std::endl;

          //(*InnerGrid)->getStorage()->recalcLeafProperty();

          this->bFirstTime = false;
        } else {
          throw generation_exception("DirichletGridConverter : buildInnerGridWithCoefs : Boundary Grid is from an unsupported grid type!");
        }
      } else {
        throw generation_exception("DirichletGridConverter : buildInnerGridWithCoefs : This method can only be called once for one instance!");
      }
    }

    void DirichletGridConverter::rebuildInnerGridWithCoefs(Grid& BoundaryGrid, DataVector& BoundaryCoefs, Grid** InnerGrid, DataVector** InnerCoefs) {
      if (this->bFirstTime == false) {
        if (strcmp(BoundaryGrid.getType(), "linearBoundary") == 0 || strcmp(BoundaryGrid.getType(), "linearTrapezoidBoundary") == 0
            || strcmp(BoundaryGrid.getType(), "linearStretchedTrapezoidBoundary") == 0) {
          GridStorage* myGridStorage = BoundaryGrid.getStorage();

          // determine the number of grid points for both grids
          this->numTotalGridPoints = myGridStorage->size();
          this->numInnerGridPoints = myGridStorage->getNumInnerPoints();

          // allocate the translation array for the coefficients
          delete[] this->conCoefArray;
          this->conCoefArray = new size_t[this->numInnerGridPoints];

          // Get the algorithmic dimensions
          std::vector<size_t> BSalgoDims = BoundaryGrid.getAlgorithmicDimensions();

          // create new inner Grid, with one grid point
          (*InnerGrid)->getStorage()->emptyStorage();

          // Set algorithmic dimensions for inner Grid
          (*InnerGrid)->setAlgorithmicDimensions(BSalgoDims);

          // create new DataVector for storing the inner grid's coefficients
          delete (*InnerCoefs);
          *InnerCoefs = new DataVector(this->numInnerGridPoints);

          // Iterate through all grid points and filter inner points
          size_t numInner = 0;

          for (size_t i = 0; i < this->numTotalGridPoints; i++) {
            GridIndex* curPoint = (*myGridStorage)[i];

            if (curPoint->isInnerPoint() == true) {
              // handle coefficients
              this->conCoefArray[numInner] = i;
              (*InnerCoefs)->set(numInner, BoundaryCoefs.get(i));
              numInner++;
              // insert point into inner grid
              (*InnerGrid)->getStorage()->insert(*curPoint);
            }
          }

          //(*InnerGrid)->getStorage()->recalcLeafProperty();
        } else {
          throw generation_exception("DirichletGridConverter : rebuildInnerGridWithCoefs : Boundary Grid is from an unsupported grid type!");
        }
      } else {
        throw generation_exception("DirichletGridConverter : rebuildInnerGridWithCoefs : This method can only be called after initial build!");
      }
    }

    void DirichletGridConverter::calcInnerCoefs(DataVector& BoundaryCoefs, DataVector& InnerCoefs) {
      if (this->numInnerGridPoints > 0) {
        for (size_t i = 0; i < this->numInnerGridPoints; i++) {
          InnerCoefs.set(i, BoundaryCoefs.get(this->conCoefArray[i]));
        }
      } else {
        throw generation_exception("DirichletGridConverter : calcInnerCoefs : The inner grid has no gridpoints! Check adaptivity calls and settings!");
      }
    }

    void DirichletGridConverter::updateBoundaryCoefs(DataVector& BoundaryCoefs, DataVector& InnerCoefs) {
      if (this->numInnerGridPoints > 0) {
        for (size_t i = 0; i < this->numInnerGridPoints; i++) {
          BoundaryCoefs.set(this->conCoefArray[i], InnerCoefs.get(i));
        }
      } else {
        throw generation_exception("DirichletGridConverter : updateBoundaryCoefs : The inner grid has no gridpoints! Check adaptivity calls and settings!");
      }
    }

  }
}
