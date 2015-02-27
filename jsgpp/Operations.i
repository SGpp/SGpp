// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

namespace sg
{
namespace base
{

class RefinementFunctor
{
public:
	typedef double value_type;

	virtual double operator()(sg::base::GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class CoarseningFunctor
{
public:
	typedef double value_type;

	virtual double operator()(sg::base::GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class GridGenerator
{
public:
	virtual void regular(size_t level) = 0;
	virtual void cliques(int level, size_t clique_size) = 0;
	virtual void truncated(size_t level,size_t l_user) = 0;
	virtual void refine(sg::base::RefinementFunctor* func) = 0;
	virtual void coarsen(sg::base::CoarseningFunctor* func, sg::base::DataVector* alpha) = 0;
	virtual void coarsenNFirstOnly(sg::base::CoarseningFunctor* func, sg::base::DataVector* alpha, size_t numFirstOnly) = 0;
	virtual size_t getNumberOfRefinablePoints() = 0;
	virtual size_t getNumberOfRemovablePoints() = 0;
	virtual void refineMaxLevel(sg::base::RefinementFunctor* func, unsigned int maxLevel) = 0;
	virtual size_t getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel) = 0;
};

class OperationMultipleEval
{
public:
	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;
	virtual void multTranspose(sg::base::DataVector& soruce, sg::base::DataVector& result) = 0;
};

class OperationMatrix
{
public:
	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;
};

class OperationConvert
{
public:
	virtual void doConvertToLinear(sg::base::DataVector& alpha) = 0;
	virtual void doConvertFromLinear(sg::base::DataVector& alpha) = 0;
};

class OperationEval
{
public:
	virtual double eval(sg::base::DataVector& alpha, sg::base::DataVector& point) = 0;
};
}

namespace datadriven {
%nodefaultdtor sg::datadriven::OperationTest;
class OperationTest
{
public:
	virtual double test(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes) = 0;
	virtual double testMSE(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& refValues) = 0;
	virtual double testWithCharacteristicNumber(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& charaNumbers) = 0;
  	virtual void calculateROCcurve(sg::base::DataVector& alpha, sg::base::DataMatrix& data, sg::base::DataVector& classes, sg::base::DataVector& thresholds, sg::base::DataMatrix& ROC_curve) = 0;
};

class OperationRegularizationDiagonal
{
public:
	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;
};

}

namespace base {
class OperationHierarchisation
{
public:
	virtual void doHierarchisation(DataVector& node_values) = 0;
	virtual void doDehierarchisation(DataVector& alpha) = 0;
};
}

namespace base {
class OperationQuadrature
{
public:
	virtual double doQuadrature(DataVector& alpha) = 0;
};
}

//-     namespace parallel --------------------------------------------
namespace parallel {

class OperationMultipleEvalVectorized
{
public:
	virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;
	virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result) = 0;
	virtual void rebuildLevelAndIndex() = 0;
};

}
//- end namespace parallel --------------------------------------------

}