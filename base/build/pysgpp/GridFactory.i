// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%{
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
%}

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

%apply std::string *OUTPUT { std::string& ostr };
%apply std::string *INPUT { std::string& istr };


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

typedef enum mail_ {
      Linear                        =  0,
      LinearStretched               =  1,
      LinearL0Boundary              =  2,
      LinearBoundary                =  3,
      LinearStretchedBoundary       =  4,
      LinearTruncatedBoundary       =  5,
      ModLinear                     =  6,
      Poly                          =  7,
      PolyBoundary                  =  8,
      ModPoly                       =  9,
      ModWavelet                    = 10,
      ModBspline                    = 11,
      Prewavelet                    = 12,
      SquareRoot                    = 13,
      Periodic                      = 14,
      LinearClenshawCurtis          = 15,
      Bspline                       = 16,
      BsplineBoundary               = 17,
      BsplineClenshawCurtis         = 18,
      Wavelet                       = 19,
      WaveletBoundary               = 20,
      FundamentalSpline             = 21,
      ModFundamentalSpline          = 22,
      ModBsplineClenshawCurtis      = 23
} GridType;

class Grid
{
public:
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
	
  static Grid* unserialize(std::string& istr);
	
protected:
  Grid();
  Grid(Grid& o);

public:
  virtual ~Grid();

public:	
  virtual SGPP::base::GridGenerator* createGridGenerator() = 0;
  virtual SGPP::base::GridStorage* getStorage();
  virtual SGPP::base::BoundingBox* getBoundingBox();
  virtual SGPP::base::Stretching* getStretching();

  virtual const char* getType() = 0;
  virtual const SBasis& getBasis() = 0;
  virtual void serialize(std::string& ostr);
  void refine(SGPP::base::DataVector* vector, int num);
  float_t eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point);
  void insertPoint(size_t dim, unsigned int levels[], unsigned int indeces[], bool isLeaf);
  int getSize();
};
}
}

//these are just two new interfaces for consistency with Memento design pattern
%extend SGPP::base::Grid{
	SGPP::base::Grid* createMemento(){
		return $self;
	}
	
	static SGPP::base::Grid* setMemento(std::string& istr){
		return SGPP::base::Grid::unserialize(istr);
	}
};

// extend the grid by a function that returns the maximum degree of the basis
// which is important for polynomials and bsplines
%extend SGPP::base::Grid{
    int getDegree() {
        if (strcmp($self->getType(), "poly") == 0) {
            return ((SGPP::base::PolyGrid*) $self)->getDegree();
        };
        if (strcmp($self->getType(), "polyBoundary") == 0) {
            return ((SGPP::base::PolyBoundaryGrid*) $self)->getDegree();
        };

        return 1;
    };
};	
