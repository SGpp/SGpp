// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "std_string.i"

%newobject SGPP::base::Grid::createLinearGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearStretchedGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject SGPP::base::Grid::createLinearClenshawCurtisGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearBoundaryGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearBoundaryGrid(SGPP::base::BoudingBox& BB);
%newobject SGPP::base::Grid::createLinearStretchedBoundaryGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearStretchedBoundaryGrid(SGPP::base::Stretching& BB);
%newobject SGPP::base::Grid::createModLinearGrid(size_t dim);
%newobject SGPP::base::Grid::createPolyGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createPolyBoundaryGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createModPolyGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createWaveletGrid(size_t dim);
%newobject SGPP::base::Grid::createWaveletBoundaryGrid(size_t dim);
%newobject SGPP::base::Grid::createModWaveletGrid(size_t dim);
%newobject SGPP::base::Grid::createBsplineGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createBsplineBoundaryGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createBsplineClenshawCurtisGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createModBsplineGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createModBsplineClenshawCurtisGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createFundamentalSplineGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createModFundamentalSlineGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createLinearTruncatedBoundaryGrid(size_t dim);
%newobject SGPP::base::Grid::createSquareRootGrid(size_t dim);
%newobject SGPP::base::Grid::createPrewaveletGrid(size_t dim);
%newobject SGPP::base::Grid::createPeriodicGrid(size_t dim);

%newobject SGPP::base::Grid::unserialize(std::string& istr);

%newobject SGPP::base::Grid::createGridGenerator();

%include "stl.i"
%include "typemaps.i"

//%apply std::string *OUTPUT { std::string& ostr };
//%apply std::string *INPUT { std::string& istr };


using namespace SGPP;
namespace SGPP
{
namespace base
{

struct RegularGridConfiguration {
      /// Grid Type, see enum
      SGPP::base::GridType type_;
      /// number of dimensions
      size_t dim_;
      /// number of levels
      int level_;
    };

struct AdpativityConfiguration {
      /// number of refinements
      size_t numRefinements_;
      /// refinement threshold for surpluses
      float_t threshold_;
      /// refinement type: false: classic, true: maxLevel
      bool maxLevelType_;
      /// max. number of points to be refined
      size_t noPoints_;
      /// max. percent of points to be refined
      float_t percent_;
    };

enum class GridType {
  Linear,                       //  0
  LinearStretched,              //  1
  LinearL0Boundary,             //  2
  LinearBoundary,               //  3
  LinearStretchedBoundary,      //  4
  LinearTruncatedBoundary,      //  5
  ModLinear,                    //  6
  Poly,                         //  7
  PolyBoundary,                 //  8
  ModPoly,                      //  9
  ModWavelet,                   // 10
  ModBspline,                   // 11
  Prewavelet,                   // 12
  SquareRoot,                   // 13
  Periodic,                     // 14
  LinearClenshawCurtis,         // 15
  Bspline,                      // 16
  BsplineBoundary,              // 17
  BsplineClenshawCurtis,        // 18
  Wavelet,                      // 19
  WaveletBoundary,              // 20
  FundamentalSpline,            // 21
  ModFundamentalSpline,         // 22
  ModBsplineClenshawCurtis,     // 23
  LinearStencil,                // 24
  ModLinearStencil              // 25
};
    
class Grid
{
protected:
  Grid();
  Grid(Grid& o);

public:
  virtual ~Grid();

public:	
  virtual SGPP::base::GridStorage& getStorage();
  virtual SGPP::base::BoundingBox& getBoundingBox();
  virtual SGPP::base::Stretching& getStretching();

  virtual SGPP::base::GridType getType() = 0;
  virtual const SBasis& getBasis() = 0;
  virtual void serialize(std::string& ostr);
  void refine(SGPP::base::DataVector* vector, int num);
  void insertPoint(size_t dim, unsigned int levels[], unsigned int indeces[], bool isLeaf);
  int getSize();
};
}
}

// SWIG doesn't support std::unique_ptr yet
// ==> use plain pointers instead (with %newobject)
%extend SGPP::base::Grid {
  static SGPP::base::Grid* createLinearGridStencil(size_t dim) {
    return SGPP::base::Grid::createLinearGridStencil(dim).release();
  }

  static SGPP::base::Grid* createModLinearGridStencil(size_t dim) {
    return SGPP::base::Grid::createModLinearGridStencil(dim).release();
  }

  static SGPP::base::Grid* createLinearGrid(size_t dim) {
    return SGPP::base::Grid::createLinearGrid(dim).release();
  }

  static SGPP::base::Grid* createLinearStretchedGrid(size_t dim) {
    return SGPP::base::Grid::createLinearStretchedGrid(dim).release();
  }

  static SGPP::base::Grid* createLinearBoundaryGrid(size_t dim, SGPP::base::level_t boundaryLevel) {
    return SGPP::base::Grid::createLinearBoundaryGrid(dim, boundaryLevel).release();
  }

  static SGPP::base::Grid* createLinearStretchedBoundaryGrid(size_t dim) {
    return SGPP::base::Grid::createLinearStretchedBoundaryGrid(dim).release();
  }

  static SGPP::base::Grid* createLinearClenshawCurtisGrid(size_t dim) {
    return SGPP::base::Grid::createLinearClenshawCurtisGrid(dim).release();
  }

  static SGPP::base::Grid* createModLinearGrid(size_t dim) {
    return SGPP::base::Grid::createModLinearGrid(dim).release();
  }

  static SGPP::base::Grid* createPolyGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createPolyGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createPolyBoundaryGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createPolyBoundaryGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createWaveletGrid(size_t dim) {
    return SGPP::base::Grid::createWaveletGrid(dim).release();
  }

  static SGPP::base::Grid* createWaveletBoundaryGrid(size_t dim) {
    return SGPP::base::Grid::createWaveletBoundaryGrid(dim).release();
  }

  static SGPP::base::Grid* createModWaveletGrid(size_t dim) {
    return SGPP::base::Grid::createModWaveletGrid(dim).release();
  }

  static SGPP::base::Grid* createBsplineGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createBsplineGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createBsplineBoundaryGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createBsplineBoundaryGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createBsplineClenshawCurtisGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createBsplineClenshawCurtisGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createModBsplineGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createModBsplineGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createModBsplineClenshawCurtisGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createModBsplineClenshawCurtisGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createFundamentalSplineGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createFundamentalSplineGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createModFundamentalSplineGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createModFundamentalSplineGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createSquareRootGrid(size_t dim) {
    return SGPP::base::Grid::createSquareRootGrid(dim).release();
  }

  static SGPP::base::Grid* createPrewaveletGrid(size_t dim) {
    return SGPP::base::Grid::createPrewaveletGrid(dim).release();
  }

  static SGPP::base::Grid* createLinearTruncatedBoundaryGrid(size_t dim) {
    return SGPP::base::Grid::createLinearTruncatedBoundaryGrid(dim).release();
  }

  static SGPP::base::Grid* createModPolyGrid(size_t dim, size_t degree) {
    return SGPP::base::Grid::createModPolyGrid(dim, degree).release();
  }

  static SGPP::base::Grid* createPeriodicGrid(size_t dim) {
    return SGPP::base::Grid::createPeriodicGrid(dim).release();
  }

  static SGPP::base::Grid* unserialize(std::string& istr) {
    return SGPP::base::Grid::unserialize(istr).release();
  }

  SGPP::base::GridGenerator* createGridGenerator() {
    return $self->createGridGenerator().release();
  }
};
  
// these are just two new interfaces for consistency with Memento design pattern
%extend SGPP::base::Grid {
  SGPP::base::Grid* createMemento() {
    return $self;
  }

  static SGPP::base::Grid* setMemento(std::string& istr) {
    return SGPP::base::Grid::unserialize(istr).release();
  }
};

/*%typemap(out) std::unique_ptr<SGPP::base::Grid> %{
    $result = BOOM;
%}*/
	
