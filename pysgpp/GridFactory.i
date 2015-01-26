/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)


%newobject SGPP::base::Grid::createLinearGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearStretchedGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearBoundaryGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearTrapezoidBoundaryGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearTrapezoidBoundaryGrid(SGPP::base::BoudingBox& BB);
%newobject SGPP::base::Grid::createLinearStretchedTrapezoidBoundaryGrid(size_t dim);
%newobject SGPP::base::Grid::createLinearStretchedTrapezoidBoundaryGrid(SGPP::base::Stretching& BB);
%newobject SGPP::base::Grid::createModLinearGrid(size_t dim);
%newobject SGPP::base::Grid::createPolyGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createModPolyGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createModWaveletGrid(size_t dim);
%newobject SGPP::base::Grid::createModBsplineGrid(size_t dim, size_t degree);
%newobject SGPP::base::Grid::createTruncatedTrapezoidGrid(size_t dim);
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
    
typedef enum mail_ {
    Linear = 0,
    LinearStretched = 1,
    LinearBoundary = 2,
    LinearTrapezoidBoundary = 3,
    LinearStretchedTrapezoidBoundary = 4,
    ModLinear = 5,
    Poly = 6,
    ModPoly = 7,
    ModWavelet = 8,
    ModBspline = 9,
    Prewavelet = 10,
    SquareRoot = 11,
    TruncatedTrapezoid = 12,
    Periodic = 13
} GridType;

class Grid
{
public:
	static Grid* createLinearGrid(size_t dim);
	static Grid* createLinearStretchedGrid(size_t dim);
	static Grid* createLinearBoundaryGrid(size_t dim);
	static Grid* createLinearTrapezoidBoundaryGrid(size_t dim);
	static Grid* createLinearStretchedTrapezoidBoundaryGrid(size_t dim);
	static Grid* createModLinearGrid(size_t dim);
	static Grid* createPolyGrid(size_t dim, size_t degree);
	static Grid* createModPolyGrid(size_t dim, size_t degree);
	static Grid* createModWaveletGrid(size_t dim);
	static Grid* createModBsplineGrid(size_t dim, size_t degree);
    static Grid* createSquareRootGrid(size_t dim);
	static Grid* createTruncatedTrapezoidGrid(size_t dim);
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
	double eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point);
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

	
