/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

%include "std_string.i"

//namespace std {
//     %typemap(in) string & = const string &;
//     // others
//}


%newobject sg::Grid::createLinearGrid(size_t dim);
%newobject sg::Grid::createLinearBoundaryGrid(size_t dim);
%newobject sg::Grid::createLinearTrapezoidBoundaryGrid(size_t dim);
%newobject sg::Grid::createLinearTrapezoidBoundaryGrid(BoudingBox& BB);
%newobject sg::Grid::createModLinearGrid(size_t dim);
%newobject sg::Grid::createPolyGrid(size_t dim, size_t degree);
%newobject sg::Grid::createModPolyGrid(size_t dim, size_t degree);
%newobject sg::Grid::createModWaveletGrid(size_t dim);
%newobject sg::Grid::createModBsplineGrid(size_t dim, size_t degree);

%newobject sg::Grid::unserialize(const std::string& istr);

%newobject sg::Grid::createOperationB();
%newobject sg::Grid::createGridGenerator();
%newobject sg::Grid::createOperationLaplace();
%newobject sg::Grid::createOperationEval();
%newobject sg::Grid::createOperationTest();
%newobject sg::Grid::createOperationHierarchisation();

%include "stl.i"
%include "typemaps.i"


//%apply std::string *OUTPUT { std::string& ostr };
//%apply std::string *INPUT { std::string& istr };



//void getMemento();


namespace sg
{

class Grid
{
public:
	static Grid* createLinearGrid(size_t dim);
	static Grid* createLinearBoundaryGrid(size_t dim);
	static Grid* createLinearTrapezoidBoundaryGrid(size_t dim);
	static Grid* createModLinearGrid(size_t dim);
	static Grid* createPolyGrid(size_t dim, size_t degree);
	static Grid* createModPolyGrid(size_t dim, size_t degree);
	static Grid* createModWaveletGrid(size_t dim);
	static Grid* createModBsplineGrid(size_t dim, size_t degree);
	
	static Grid* unserialize(const std::string& istr);
	
protected:
	Grid();
	Grid(Grid& o);

public:
	virtual ~Grid();

public:	
	virtual GridGenerator* createGridGenerator() = 0;
	virtual OperationB* createOperationB() = 0;
	virtual OperationEval* createOperationEval() = 0;
	virtual OperationTest* createOperationTest() = 0;
	virtual OperationMatrix* createOperationLaplace() = 0;
	virtual OperationMatrix* createOperationIdentity() = 0;
	virtual OperationHierarchisation* createOperationHierarchisation() = 0;
	
	// @todo remove this when done
	virtual OperationMatrix* createOperationUpDownTest() = 0;
	
	virtual OperationMatrix* createOperationDelta(DataVector& coef) = 0;
	virtual OperationMatrix* createOperationGamma(DataVector& coef) = 0;
	
	virtual GridStorage* getStorage();
	virtual BoundingBox* getBoundingBox();
	virtual void setBoundingBox(BoundingBox& bb);

	virtual const char* getType() = 0;	
	virtual void serialize(std::string& ostr);
	virtual std::string serialize();
	void refine(DataVector* vector, int num);
	double eval(DataVector& alpha, DataVector& point);
	void insertPoint(size_t dim, unsigned int levels[], unsigned int indeces[], bool isLeaf);
	int getSize();
	
};



}

//these are just two new interfaces for consistency with Memento design pattern
%extend sg::Grid{
	Grid* createMemento(){
		return $self;
	}
	
	static sg::Grid* setMemento(std::string& istr){
		return sg::Grid::unserialize(istr);
	}
};

	
