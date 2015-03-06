// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

namespace SGPP
{
namespace base
{

class RefinementFunctor
{
public:
    typedef float_t value_type;

    virtual float_t operator()(SGPP::base::GridStorage* storage, size_t seq) = 0;
    virtual float_t start() = 0;    
};

class CoarseningFunctor
{
public:
    typedef float_t value_type;

    virtual float_t operator()(SGPP::base::GridStorage* storage, size_t seq) = 0;
    virtual float_t start() = 0;    
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
    virtual float_t eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point) = 0;
};

class OperationNaiveEval
{
public:
    virtual float_t eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point) = 0;
};

class OperationNaiveEvalGradient
{
public:
    virtual float_t eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point, SGPP::base::DataVector& gradient) = 0;
};

class OperationNaiveEvalHessian
{
public:
    virtual float_t eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point, SGPP::base::DataVector& gradient, SGPP::base::DataMatrix& hessian) = 0;
};

class OperationNaiveEvalPartialDerivative
{
public:
    virtual float_t eval(SGPP::base::DataVector& alpha, SGPP::base::DataVector& point, size_t deriv_dim) = 0;
};

class OperationHierarchisation
{
public:
    virtual void doHierarchisation(SGPP::base::DataVector& node_values) = 0;
    virtual void doDehierarchisation(SGPP::base::DataVector& alpha) = 0;
};

class OperationQuadrature
{
public:
    virtual float_t doQuadrature(SGPP::base::DataVector& alpha) = 0;
};

}

//-     namespace parallel --------------------------------------------
/*namespace parallel {

class OperationMultipleEvalVectorized
{
public:
	virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;
	virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result) = 0;
	virtual void rebuildLevelAndIndex() = 0;
};

}*/
//- end namespace parallel --------------------------------------------

}