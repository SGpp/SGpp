// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "std_string.i"

%newobject sgpp::base::Grid::createGrid(RegularGridConfiguration gridConfig);
%newobject sgpp::base::Grid::createLinearGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearStretchedGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearBoundaryGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject sgpp::base::Grid::createLinearClenshawCurtisGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearClenshawCurtisBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject sgpp::base::Grid::createModLinearClenshawCurtisGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearBoundaryGrid(sgpp::base::BoudingBox& BB);
%newobject sgpp::base::Grid::createLinearStretchedBoundaryGrid(size_t dim);
%newobject sgpp::base::Grid::createLinearStretchedBoundaryGrid(sgpp::base::Stretching& BB);
%newobject sgpp::base::Grid::createModLinearGrid(size_t dim);
%newobject sgpp::base::Grid::createPolyGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createPolyBoundaryGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createPolyBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel);
%newobject sgpp::base::Grid::createModPolyGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createWaveletGrid(size_t dim);
%newobject sgpp::base::Grid::createWaveletBoundaryGrid(size_t dim);
%newobject sgpp::base::Grid::createWaveletBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject sgpp::base::Grid::createModWaveletGrid(size_t dim);
%newobject sgpp::base::Grid::createBsplineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createBsplineBoundaryGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createBsplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel);
%newobject sgpp::base::Grid::createBsplineClenshawCurtisGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createBsplineClenshawCurtisGrid(size_t dim, size_t degree, size_t boundaryLevel);
%newobject sgpp::base::Grid::createModBsplineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModBsplineClenshawCurtisGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createFundamentalSplineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModFundamentalSlineGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createLinearTruncatedBoundaryGrid(size_t dim);
%newobject sgpp::base::Grid::createSquareRootGrid(size_t dim);
%newobject sgpp::base::Grid::createPrewaveletGrid(size_t dim);
%newobject sgpp::base::Grid::createPeriodicGrid(size_t dim);
%newobject sgpp::base::Grid::createPolyClenshawCurtisBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel);
%newobject sgpp::base::Grid::createPolyClenshawCurtisGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createModPolyClenshawCurtisGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(size_t dim, size_t degree);
%newobject sgpp::base::Grid::createNaturalBsplineBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject sgpp::base::Grid::createNakBsplineBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject sgpp::base::Grid::createModNakBsplineGrid(size_t dim);
%newobject sgpp::base::Grid::createWeaklyFundamentalSplineBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject sgpp::base::Grid::createWeaklyFundamentalNakSplineBoundaryGrid(size_t dim, size_t boundaryLevel);
%newobject sgpp::base::Grid::createModWeaklyFundamentalNakSplineGrid(size_t dim);
%newobject sgpp::base::Grid::createFundamentalSplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel);
%newobject sgpp::base::Grid::createFundamentalNakSplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel);

%newobject sgpp::base::Grid::unserialize(std::string& istr);
%newobject sgpp::base::Grid::createGridOfEquivalentType(size_t numDims);
%newobject sgpp::base::Grid::clone();

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
  /// max. polynomial degree for poly basis
  size_t maxDegree_;
  /// level of boundary grid
  sgpp::base::level_t boundaryLevel_;
  /// string to serialized grid
  std::string filename_;
    };

struct AdaptivityConfiguration {
    /// number of refinements
    size_t numRefinements_;
    /// refinement threshold
    double refinementThreshold_;
    /// coarsening threshold
    double coarseningThreshold_;
    /// refinement type: false: classic, true: maxLevel
    bool maxLevelType_;
    /// max. number of points to be refined
    size_t numRefinementPoints_;
    /// max. number of points to be coarsened
    size_t numCoarseningPoints_;
    /// max. percent of points to be refined
    double percent_;
    /// other refinement strategy, that is more expensive, but yields better results
    bool errorBasedRefinement_ = false;
    /// prevent coarsening of initial grid points, needed for some decompositions
    bool coarsenInitialPoints_;
    };

enum class GridType {
  Linear,                                   //  0
  LinearStretched,                          //  1
  LinearL0Boundary,                         //  2
  LinearBoundary,                           //  3
  LinearStretchedBoundary,                  //  4
  LinearTruncatedBoundary,                  //  5
  ModLinear,                                //  6
  Poly,                                     //  7
  PolyBoundary,                             //  8
  ModPoly,                                  //  9
  ModWavelet,                               // 10
  ModBspline,                               // 11
  Prewavelet,                               // 12
  SquareRoot,                               // 13
  Periodic,                                 // 14
  LinearClenshawCurtisBoundary,             // 15
  Bspline,                                  // 16
  BsplineBoundary,                          // 17
  BsplineClenshawCurtis,                    // 18
  Wavelet,                                  // 19
  WaveletBoundary,                          // 20
  FundamentalSpline,                        // 21
  ModFundamentalSpline,                     // 22
  ModBsplineClenshawCurtis,                 // 23
  LinearStencil,                            // 24
  ModLinearStencil,                         // 25
  PolyClenshawCurtisBoundary,               // 26
  PolyClenshawCurtis,                       // 27
  LinearClenshawCurtis,                     // 28
  ModPolyClenshawCurtis,                    // 29
  ModLinearClenshawCurtis,                  // 30
  NakBsplineBoundaryCombigrid,              // 31
  NaturalBsplineBoundary,                   // 32
  NakBsplineBoundary,                  // 33
  ModNakBspline,                       // 34
  WeaklyFundamentalSplineBoundary,          // 35
  WeaklyFundamentalNakSplineBoundary,  // 36
  ModWeaklyFundamentalNakSpline,       // 37
  FundamentalSplineBoundary,                // 38
  FundamentalNakSplineBoundary,        // 39
};

class Grid
{
public:
  static Grid* createGrid(sgpp::base::RegularGridConfiguration gridConfig);
  static Grid* createLinearGrid(size_t dim);
  static Grid* createLinearStretchedGrid(size_t dim);
  static Grid* createLinearBoundaryGrid(size_t dim, size_t boundaryLevel=1);
  static Grid* createLinearClenshawCurtisGrid(size_t dim);
  static Grid* createLinearClenshawCurtisBoundaryGrid(size_t dim, size_t boundaryLevel=1);
  static Grid* createModLinearClenshawCurtisGrid(size_t dim);
  static Grid* createLinearStretchedBoundaryGrid(size_t dim);
  static Grid* createModLinearGrid(size_t dim);
  static Grid* createPolyGrid(size_t dim, size_t degree);
  static Grid* createPolyBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
  static Grid* createModPolyGrid(size_t dim, size_t degree);
  static Grid* createWaveletGrid(size_t dim);
  static Grid* createWaveletBoundaryGrid(size_t dim, size_t boundaryLevel=1);
  static Grid* createModWaveletGrid(size_t dim);
  static Grid* createBsplineGrid(size_t dim, size_t degree);
  static Grid* createBsplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
  static Grid* createBsplineClenshawCurtisGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
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
  static Grid* createPolyClenshawCurtisBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
  static Grid* createPolyClenshawCurtisGrid(size_t dim, size_t degree);
  static Grid* createModPolyClenshawCurtisGrid(size_t dim, size_t degree);
  static Grid* createNakBsplineBoundaryCombigridGrid(size_t dim, size_t degree);
  static Grid* createNaturalBsplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
  static Grid* createNakBsplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
  static Grid* createModNakBsplineGrid(size_t dim, size_t degree);
  static Grid* createWeaklyFundamentalSplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
  static Grid* createWeaklyFundamentalNakSplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
  static Grid* createFundamentalSplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);
  static Grid* createFundamentalNakSplineBoundaryGrid(size_t dim, size_t degree, size_t boundaryLevel=1);

  static Grid* unserialize(const std::string& istr);

  static sgpp::base::GridType stringToGridType(const std::string& gridType);
	
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
  
  std::string getTypeAsString();

  Grid* createGridOfEquivalentType(size_t numDims);
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
    size_t getDegree() {
	return $self->getBasis().getDegree();
    };
};
