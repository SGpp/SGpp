// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "std_string.i"

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

%newobject sgpp::base::Grid::unserialize(std::string& istr);

%include "stl.i"
%include "typemaps.i"

//%apply std::string *OUTPUT { std::string& ostr };
//%apply std::string *INPUT { std::string& istr };


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
};
}
}

// SWIG doesn't support std::unique_ptr yet
// ==> use plain pointers instead (with %newobject)
%extend sgpp::base::Grid {
  static sgpp::base::Grid* createLinearGridStencil(size_t dim) {
    return sgpp::base::Grid::createLinearGridStencil(dim).release();
  }

  static sgpp::base::Grid* createModLinearGridStencil(size_t dim) {
    return sgpp::base::Grid::createModLinearGridStencil(dim).release();
  }

  static sgpp::base::Grid* createLinearGrid(size_t dim) {
    return sgpp::base::Grid::createLinearGrid(dim).release();
  }

  static sgpp::base::Grid* createLinearStretchedGrid(size_t dim) {
    return sgpp::base::Grid::createLinearStretchedGrid(dim).release();
  }

  static sgpp::base::Grid* createLinearBoundaryGrid(size_t dim, sgpp::base::level_t boundaryLevel) {
    return sgpp::base::Grid::createLinearBoundaryGrid(dim, boundaryLevel).release();
  }

  static sgpp::base::Grid* createLinearStretchedBoundaryGrid(size_t dim) {
    return sgpp::base::Grid::createLinearStretchedBoundaryGrid(dim).release();
  }

  static sgpp::base::Grid* createLinearClenshawCurtisGrid(size_t dim) {
    return sgpp::base::Grid::createLinearClenshawCurtisGrid(dim).release();
  }

  static sgpp::base::Grid* createModLinearGrid(size_t dim) {
    return sgpp::base::Grid::createModLinearGrid(dim).release();
  }

  static sgpp::base::Grid* createPolyGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createPolyGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createPolyBoundaryGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createPolyBoundaryGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createWaveletGrid(size_t dim) {
    return sgpp::base::Grid::createWaveletGrid(dim).release();
  }

  static sgpp::base::Grid* createWaveletBoundaryGrid(size_t dim) {
    return sgpp::base::Grid::createWaveletBoundaryGrid(dim).release();
  }

  static sgpp::base::Grid* createModWaveletGrid(size_t dim) {
    return sgpp::base::Grid::createModWaveletGrid(dim).release();
  }

  static sgpp::base::Grid* createBsplineGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createBsplineGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createBsplineBoundaryGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createBsplineBoundaryGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createBsplineClenshawCurtisGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createBsplineClenshawCurtisGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createModBsplineGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createModBsplineGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createModBsplineClenshawCurtisGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createModBsplineClenshawCurtisGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createFundamentalSplineGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createFundamentalSplineGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createModFundamentalSplineGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createModFundamentalSplineGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createSquareRootGrid(size_t dim) {
    return sgpp::base::Grid::createSquareRootGrid(dim).release();
  }

  static sgpp::base::Grid* createPrewaveletGrid(size_t dim) {
    return sgpp::base::Grid::createPrewaveletGrid(dim).release();
  }

  static sgpp::base::Grid* createLinearTruncatedBoundaryGrid(size_t dim) {
    return sgpp::base::Grid::createLinearTruncatedBoundaryGrid(dim).release();
  }

  static sgpp::base::Grid* createModPolyGrid(size_t dim, size_t degree) {
    return sgpp::base::Grid::createModPolyGrid(dim, degree).release();
  }

  static sgpp::base::Grid* createPeriodicGrid(size_t dim) {
    return sgpp::base::Grid::createPeriodicGrid(dim).release();
  }

  static sgpp::base::Grid* unserialize(std::string& istr) {
    return sgpp::base::Grid::unserialize(istr).release();
  }
};
  
// these are just two new interfaces for consistency with Memento design pattern
%extend sgpp::base::Grid {
  sgpp::base::Grid* createMemento() {
    return $self;
  }

  static sgpp::base::Grid* setMemento(std::string& istr) {
    return sgpp::base::Grid::unserialize(istr).release();
  }
};

