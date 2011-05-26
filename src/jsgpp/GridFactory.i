/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

%include "std_string.i"

%newobject sg::base::Grid::createLinearGrid(size_t dim);
%newobject sg::base::Grid::createLinearStretchedGrid(size_t dim);
%newobject sg::base::Grid::createLinearBoundaryGrid(size_t dim);
%newobject sg::base::Grid::createLinearTrapezoidBoundaryGrid(size_t dim);
%newobject sg::base::Grid::createLinearTrapezoidBoundaryGrid(BoudingBox& BB);
%newobject sg::base::Grid::createLinearStretchedTrapezoidBoundaryGrid(size_t dim);
%newobject sg::base::Grid::createLinearStretchedTrapezoidBoundaryGrid(Stretching& BB);
%newobject sg::base::Grid::createModLinearGrid(size_t dim);
%newobject sg::base::Grid::createPolyGrid(size_t dim, size_t degree);
%newobject sg::base::Grid::createModPolyGrid(size_t dim, size_t degree);
%newobject sg::base::Grid::createModWaveletGrid(size_t dim);
%newobject sg::base::Grid::createModBsplineGrid(size_t dim, size_t degree);
%newobject sg::base::Grid::createTruncatedTrapezoidGrid(size_t dim);
%newobject sg::base::Grid::createSquareRootGrid(size_t dim);
%newobject sg::base::Grid::createPrewaveletGrid(size_t dim);

%newobject sg::base::Grid::unserialize(std::string& istr);

%newobject sg::base::Grid::createGridGenerator();

%include "stl.i"
%include "typemaps.i"

%apply std::string *OUTPUT { std::string& ostr };
%apply std::string *INPUT { std::string& istr };



//void getMemento();

using namespace sg;
namespace sg
{
namespace base
{

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
	
	static Grid* unserialize(const std::string& istr); 
	
protected:
	Grid();
	Grid(Grid& o);

public:
	virtual ~Grid();

public:	
	virtual GridGenerator* createGridGenerator() = 0;
	//virtual OperationMultipleEval* createOperationMultipleEval(DataMatrix* dataset) = 0;
	//virtual OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(const std::string& VecType, DataMatrix* dataset) = 0;
	//virtual OperationEval* createOperationEval() = 0;
	//virtual sg::datadriven::OperationTest* createOperationTest() = 0;
	//virtual OperationMatrix* createOperationLaplace() = 0;
	//virtual OperationMatrix* createOperationLTwoDotProduct() = 0;
	//virtual OperationMatrix* createOperationIdentity() = 0;
	//virtual OperationHierarchisation* createOperationHierarchisation() = 0;
	//virtual OperationConvert* createOperationConvert() = 0;
	
	//virtual OperationMatrix* createOperationDelta(DataVector& coef) = 0;
	//virtual OperationMatrix* createOperationGamma(DataMatrix& coef) = 0;
	
	virtual GridStorage* getStorage();
	virtual BoundingBox* getBoundingBox();
	virtual Stretching* getStretching();

	virtual const char* getType() = 0;	
	virtual void serialize(std::string& ostr); 
	void refine(DataVector* vector, int num);
	double eval(DataVector& alpha, DataVector& point);
	void insertPoint(size_t dim, unsigned int levels[], unsigned int indeces[], bool isLeaf);
	int getSize();
	
};
}
}

//these are just two new interfaces for consistency with Memento design pattern
%extend sg::base::Grid{
	Grid* createMemento(){
		return $self;
	}
	
	static Grid* setMemento(std::string& istr){
		return Grid::unserialize(istr);
	}
};

	
