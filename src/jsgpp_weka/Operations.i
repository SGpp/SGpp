// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

namespace sg
{

class RefinementFunctor
{
public:
	typedef double value_type;

	virtual double operator()(GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class CoarseningFunctor
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
	virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) = 0;
	virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly) = 0;
	virtual int getNumberOfRefinablePoints() = 0;
	virtual int getNumberOfRemovablePoints() = 0;
	virtual void refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel) = 0;
	virtual int getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel) = 0;
};

class OperationMultipleEval
{
public:
	virtual void mult(DataVector& alpha, DataVector& result) = 0;
	virtual void multTranspose(DataVector& soruce, DataVector& result) = 0;
};

class OperationMultipleEvalVectorized
{
public:
	virtual double multVectorized(DataVector& alpha, DataVector& result) = 0;
	virtual double multTransposeVectorized(DataVector& source, DataVector& result) = 0;
	virtual void rebuildLevelAndIndex() = 0;
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
	virtual double testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues) = 0;
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers) = 0;
};

class OperationHierarchisation
{
public:
	virtual void doHierarchisation(DataVector& node_values) = 0;
	virtual void doDehierarchisation(DataVector& alpha) = 0;
};

}