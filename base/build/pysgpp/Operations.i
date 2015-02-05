/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)


namespace SGPP
{
//-     namespace base ------------------------------------------------
namespace base
{

class RefinementFunctor
{
public:
	typedef double value_type;

	virtual double operator()(SGPP::base::GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class CoarseningFunctor
{
public:
	typedef double value_type;

	virtual double operator()(SGPP::base::GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class GridGenerator
{
public:
	virtual void regular(size_t level) = 0;
	virtual void cliques(int level, size_t clique_size) = 0;
	virtual void full(size_t level) = 0;
	virtual void truncated(size_t level,size_t l_user) = 0;
	virtual void refine(SGPP::base::RefinementFunctor* func) = 0;
	virtual void coarsen(SGPP::base::CoarseningFunctor* func, SGPP::base::DataVector* alpha) = 0;
	virtual void coarsenNFirstOnly(SGPP::base::CoarseningFunctor* func, SGPP::base::DataVector* alpha, size_t numFirstOnly) = 0;
	virtual int getNumberOfRefinablePoints() = 0;
	virtual int getNumberOfRemovablePoints() = 0;
	virtual void refineMaxLevel(SGPP::base::RefinementFunctor* func, unsigned int maxLevel) = 0;
	virtual int getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel) = 0;
};

class OperationMultipleEval
{
public:
	virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;
	virtual void multTranspose(SGPP::base::DataVector& soruce, SGPP::base::DataVector& result) = 0;
};

class OperationMatrix
{
public:
	virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;
};*/
/*
class OperationConvert
{
public:
	virtual void doConvertToLinear(SGPP::base::DataVector& alpha) = 0;
	virtual void doConvertFromLinear(SGPP::base::DataVector& alpha) = 0;
};

class OperationEval
{
public:
	virtual double eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point) = 0;
};

class OperationNaiveEval
{
public:
    virtual double eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point) = 0;
};

class OperationNaiveEvalGradient
{
public:
    virtual double eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point, SGPP::base::DataVector& gradient) = 0;
};

class OperationNaiveEvalHessian
{
public:
    virtual double eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point, SGPP::base::DataVector& gradient, SGPP::base::DataMatrix& hessian) = 0;
};

class OperationNaiveEvalPartialDerivative
{
public:
    virtual double eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point, size_t deriv_dim) = 0;
};

class OperationHierarchisation
{
public:
	virtual void doHierarchisation(SGPP::base::DataVector& node_values) = 0;
	virtual void doDehierarchisation(SGPP::base::DataVector& alpha) = 0;
};


}

namespace base {
class OperationQuadrature
{
public:
	virtual double doQuadrature(SGPP::base::DataVector& alpha) = 0;
};

}
//- end namespace base ------------------------------------------------

}
