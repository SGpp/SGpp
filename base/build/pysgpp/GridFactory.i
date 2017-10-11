// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%{
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
%}

%newobject sgpp::base::Grid::createGrid(RegularGridConfiguratio gridConfig);
%newobject sgpp::base::Grid::createLinearGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearStretchedGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject sgpp::base::Grid::createLinearClenshawCurtisGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearBoundaryGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearBoundaryGrid(sgpp::base::BoudingBox& BB);
%newobject sgpp::base::Grid::createLinearStretchedBoundaryGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearStretchedBoundaryGrid(sgpp::base::Stretching& BB);
%newobject sgpp::base::Grid::createModLinearGrid(size_t dim);
%newobject sgpp::base::Grid::createPolyGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createPolyBoundaryGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModPolyGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createWaveletGrid(size_t dim);
%newobject sgpp::base::Grid::createWaveletBoundaryGrid(size_t dim);
%newobject sgpp::base::Grid::createModWaveletGrid(size_t dim);
%newobject sgpp::base::Grid::createBsplineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createBsplineBoundaryGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createBsplineClenshawCurtisGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModBsplineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModBsplineClenshawCurtisGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createFundamentalSplineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModFundamentalSlineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createLinearTruncatedBoundaryGrid(size_t dim);
%newobject sgpp::base::Grid::createSquareRootGrid(size_t dim);
%newobject sgpp::base::Grid::createPrewaveletGrid(size_t dim);
%newobject sgpp::base::Grid::createPeriodicGrid(size_t dim);
%newobject sgpp::base::Grid::createNaturalBsplineBoundaryGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createNotAKnotBsplineBoundaryGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModNotAKnotBsplineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createLagrangeSplineBoundaryGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createLagrangeNotAKnotSplineBoundaryGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModLagrangeNotAKnotSplineGrid(size_t dim, size_t degree);

%newobject sgpp::base::Grid::unserialize(std::string& istr);
%newobject sgpp::base::Grid::clone();

%include "stl.i"
%include "typemaps.i"

%apply std::string *OUTPUT { std::string& ostr };
%apply std::string *INPUT { std::string& istr };


using namespace sgpp;
namespace sgpp
{
namespace base
{
struct RegularGridConfiguration {
      /// Grid Type, see enum
      sgpp::base::GridType type_;
      /// number of dimensions
      size_t dim_;
      /// number of levels
      int level_;
      /// max. polynomial degree for poly basis
      size_t maxDegree_;
      /// level of boundary grid
      size_t boundaryLevel_;
      /// subgrid selection value t
      double t_ = 0.0;
    };


struct AdpativityConfiguration {
      /// number of refinements
      size_t numRefinements_;
      /// refinement threshold for surpluses
      double threshold_;
      /// refinement type: false: classic, true: maxLevel
      bool maxLevelType_;
      /// max. number of points to be refined
      size_t noPoints_;
      /// max. percent of points to be refined
      double percent_;
    };

enum class GridType {
  Linear,                           //  0
  LinearStretched,                  //  1
  LinearL0Boundary,                 //  2
  LinearBoundary,                   //  3
  LinearStretchedBoundary,          //  4
  LinearTruncatedBoundary,          //  5
  ModLinear,                        //  6
  Poly,                             //  7
  PolyBoundary,                     //  8
  ModPoly,                          //  9
  ModWavelet,                       // 10
  ModBspline,                       // 11
  Prewavelet,                       // 12
  SquareRoot,                       // 13
  Periodic,                         // 14
  LinearClenshawCurtis,             // 15
  Bspline,                          // 16
  BsplineBoundary,                  // 17
  BsplineClenshawCurtis,            // 18
  Wavelet,                          // 19
  WaveletBoundary,                  // 20
  FundamentalSpline,                // 21
  ModFundamentalSpline,             // 22
  ModBsplineClenshawCurtis,         // 23
  LinearStencil,                    // 24
  ModLinearStencil,                 // 25
  NaturalBsplineBoundary,           // 26
  NotAKnotBsplineBoundary,          // 27
  ModNotAKnotBspline,               // 28
  LagrangeSplineBoundary,           // 29
  LagrangeNotAKnotSplineBoundary,   // 30
  ModLagrangeNotAKnotSpline,        // 31
};

class Grid
{
public:
  static Grid* createGrid(RegularGridConfiguration gridConfig);
  static Grid* createLinearGrid(size_t dim);
  static Grid* createLinearStretchedGrid(size_t dim);
  static Grid* createLinearBoundaryGrid(size_t dim, size_t boundaryLevel);
  static Grid* createLinearClenshawCurtisGrid(size_t dim);
  static Grid* createLinearBoundaryGrid(size_t dim);
  static Grid* createLinearStretchedBoundaryGrid(size_t dim);
  static Grid* createModLinearGrid(size_t dim);
  static Grid* createPolyGrid(size_t dim, size_t degree);
  static Grid* createPolyBoundaryGrid(size_t dim, size_t degree);
  static Grid* createModPolyGrid(size_t dim, size_t degree);
  static Grid* createWaveletGrid(size_t dim);
  static Grid* createWaveletBoundaryGrid(size_t dim);
  static Grid* createModWaveletGrid(size_t dim);
  static Grid* createBsplineGrid(size_t dim, size_t degree);
  static Grid* createBsplineBoundaryGrid(size_t dim, size_t degree);
  static Grid* createBsplineClenshawCurtisGrid(size_t dim, size_t degree);
  static Grid* createModBsplineGrid(size_t dim, size_t degree);
  static Grid* createModBsplineClenshawCurtisGrid(size_t dim, size_t degree);
  static Grid* createFundamentalSplineGrid(size_t dim, size_t degree);
  static Grid* createModFundamentalSplineGrid(size_t dim, size_t degree);
  static Grid* createSquareRootGrid(size_t dim);
  static Grid* createLinearTruncatedBoundaryGrid(size_t dim);
  static Grid* createPrewaveletGrid(size_t dim);
  static Grid* createLinearGridStencil(size_t dim);
  static Grid* createModLinearGridStencil(size_t dim);
  static Grid* createPeriodicGrid(size_t dim);
  static Grid* createNaturalBsplineBoundaryGrid(size_t dim, size_t degree);
  static Grid* createNotAKnotBsplineBoundaryGrid(size_t dim, size_t degree);
  static Grid* createModNotAKnotBsplineGrid(size_t dim, size_t degree);
  static Grid* createLagrangeSplineBoundaryGrid(size_t dim, size_t degree);
  static Grid* createLagrangeNotAKnotSplineBoundaryGrid(size_t dim, size_t degree);
  static Grid* createModLagrangeNotAKnotSplineGrid(size_t dim, size_t degree);
	
  static Grid* unserialize(std::string& istr);
	
protected:
  Grid();
  Grid(Grid& o);

public:
  virtual ~Grid();

public:	
  virtual sgpp::base::GridStorage& getStorage();
  virtual sgpp::base::BoundingBox& getBoundingBox();
  virtual sgpp::base::Stretching& getStretching();

  virtual sgpp::base::GridGenerator& getGenerator() = 0;
  virtual sgpp::base::GridType getType() = 0;
  virtual const SBasis& getBasis() = 0;
  virtual void serialize(std::string& ostr);
  void refine(sgpp::base::DataVector& vector, int num);
  void insertPoint(size_t dim, unsigned int levels[], unsigned int indeces[], bool isLeaf);
  int getSize();
  
  Grid* clone();
};
}
}

// these are just two new interfaces for consistency with Memento design pattern
%extend sgpp::base::Grid {
  sgpp::base::Grid* createMemento() {
    return $self;
  }

  static sgpp::base::Grid* setMemento(std::string& istr) {
    return sgpp::base::Grid::unserialize(istr);
  }
};

// extend the grid by a function that returns the maximum degree of the basis
// which is important for polynomials and bsplines
%extend sgpp::base::Grid{
    int getDegree() {
        if ($self->getType() == sgpp::base::GridType::Poly) {
            return ((sgpp::base::PolyGrid*) $self)->getDegree();
        };
        if ($self->getType() == sgpp::base::GridType::PolyBoundary) {
            return ((sgpp::base::PolyBoundaryGrid*) $self)->getDegree();
        };
        if ($self->getType() == sgpp::base::GridType::ModPoly) {
            return ((sgpp::base::ModPolyGrid*) $self)->getDegree();
        };
        return 1;
    };
};
