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
	virtual void multTranspose(SGPP::base::DataVector& source, SGPP::base::DataVector& result) = 0;
};

class OperationMatrix
{
public:
	virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;
};

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

#ifdef SG_DATADRIVEN
//-     namespace datadriven ------------------------------------------
namespace datadriven {
%nodefaultdtor SGPP::datadriven::OperationTest;
class OperationTest
{
public:
  virtual double test(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes) = 0;
	virtual double testMSE(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& refValues) = 0;
	virtual double testWithCharacteristicNumber(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& charaNumbers) = 0;
	virtual void calculateROCcurve(SGPP::base::DataVector& alpha, SGPP::base::DataMatrix& data, SGPP::base::DataVector& classes, SGPP::base::DataVector& thresholds, SGPP::base::DataMatrix& ROC_curve) = 0;
};

class OperationRegularizationDiagonal
{
public:
	virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;
	static const int HKMIX = 1;
	static const int H0HKLAPLACE = 2;
	static const int ISOTROPIC_PENALTY = 3;
	static const int ANISOTROPIC_PENALTY = 4;
};

class OperationRosenblattTransformation
{
public:
	virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* points, base::DataMatrix* pointscdf) = 0;
	virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* points, base::DataMatrix* pointscdf, size_t dim_start) = 0;
};

class OperationInverseRosenblattTransformation
{
public:
	virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* pointscdf, base::DataMatrix* points) = 0;
	virtual void doTransformation(base::DataVector* alpha, base::DataMatrix* pointscdf, base::DataMatrix* points, size_t dim_start) = 0;
};

}
//- end namespace datadriven ------------------------------------------
#endif

#ifdef SG_PARALLEL
//-     namespace parallel --------------------------------------------
namespace parallel {

class OperationMultipleEvalVectorized
{
public:
	virtual double multVectorized(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) = 0;
	virtual double multTransposeVectorized(SGPP::base::DataVector& source, SGPP::base::DataVector& result) = 0;
	virtual void rebuildLevelAndIndex() = 0;
};

}
//- end namespace parallel --------------------------------------------
#endif
}

