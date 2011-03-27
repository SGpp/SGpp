/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

namespace sg
{

class RefinementFunctor
{
public:
	typedef double value_type;

	virtual double operator()(GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class GridGenerator
{
public:
	virtual void regular(size_t level) = 0;
	virtual void refine(RefinementFunctor* func) = 0;
	virtual int getNumberOfRefinablePoints() = 0;
};

class OperationB
{
public:
	virtual void mult(DataVector& alpha, DataMatrix& data, DataVector& result) = 0;
	virtual void multTranspose(DataVector& alpha, DataMatrix& data, DataVector& result) = 0;
};

class OperationMatrix
{
public:
	virtual void mult(DataVector& alpha, DataVector& result) = 0;
};

class OperationEval
{
public:
	virtual double eval(DataVector& alpha, DataVector& point) = 0;
};

class OperationTest
{
public:
	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes) = 0;
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers) = 0;
};

class OperationHierarchisation
{
public:
	virtual void doHierarchisation(DataVector& node_values) = 0;
	virtual void doDehierarchisation(DataVector& alpha) = 0;
};

}
