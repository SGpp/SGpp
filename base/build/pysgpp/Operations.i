// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


namespace sgpp
{
//-     namespace base ------------------------------------------------
namespace base
{

class RefinementFunctor
{
public:
	typedef double value_type;

	virtual double operator()(sgpp::base::GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class CoarseningFunctor
{
public:
	typedef double value_type;

	virtual double operator()(sgpp::base::GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class GridGenerator
{
public:
	virtual void regular(size_t level) = 0;
	virtual void cliques(int level, size_t clique_size) = 0;
	virtual void full(size_t level) = 0;
	virtual void truncated(size_t level,size_t l_user) = 0;
	virtual void refine(sgpp::base::RefinementFunctor* func) = 0;
	virtual void coarsen(sgpp::base::CoarseningFunctor* func, sgpp::base::DataVector* alpha) = 0;
	virtual void coarsenNFirstOnly(sgpp::base::CoarseningFunctor* func, sgpp::base::DataVector* alpha, size_t numFirstOnly) = 0;
	virtual int getNumberOfRefinablePoints() = 0;
	virtual int getNumberOfRemovablePoints() = 0;
	virtual void refineMaxLevel(sgpp::base::RefinementFunctor* func, unsigned int maxLevel) = 0;
	virtual int getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel) = 0;
};

class OperationMultipleEval
{
public:
	virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) = 0;
	virtual void multTranspose(sgpp::base::DataVector& soruce, sgpp::base::DataVector& result) = 0;
};

class OperationMatrix
{
public:
	virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) = 0;
};*/
/*
class OperationConvert
{
public:
	virtual void doConvertToLinear(sgpp::base::DataVector& alpha) = 0;
	virtual void doConvertFromLinear(sgpp::base::DataVector& alpha) = 0;
};

class OperationEval
{
public:
	virtual double eval(sgpp::base::DataVector& alpha, sgpp::base::DataVector& point) = 0;
};

class OperationNaiveEval
{
public:
    virtual double eval(sgpp::base::DataVector& alpha, sgpp::base::DataVector& point) = 0;
};

class OperationNaiveEvalGradient
{
public:
    virtual double eval(sgpp::base::DataVector& alpha, sgpp::base::DataVector& point, sgpp::base::DataVector& gradient) = 0;
};

class OperationNaiveEvalHessian
{
public:
    virtual double eval(sgpp::base::DataVector& alpha, sgpp::base::DataVector& point, sgpp::base::DataVector& gradient, sgpp::base::DataMatrix& hessian) = 0;
};

class OperationNaiveEvalPartialDerivative
{
public:
    virtual double eval(sgpp::base::DataVector& alpha, sgpp::base::DataVector& point, size_t deriv_dim) = 0;
};

class OperationHierarchisation
{
public:
	virtual void doHierarchisation(sgpp::base::DataVector& node_values) = 0;
	virtual void doDehierarchisation(sgpp::base::DataVector& alpha) = 0;
};

class OperationQuadrature
{
public:
	virtual double doQuadrature(sgpp::base::DataVector& alpha) = 0;
};

}
//- end namespace base ------------------------------------------------

}
