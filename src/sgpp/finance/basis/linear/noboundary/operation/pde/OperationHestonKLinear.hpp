/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#ifndef OPERATIONHESTONKLINEAR_HPP
#define OPERATIONHESTONKLINEAR_HPP

#include "pde/algorithm/UpDownFourOpDims.hpp"

namespace sg
{
namespace finance
{

/**
 * Implements the Heston B-Operation (corresponds to matrix B in Master's thesis), that is needed
 * the solve the multidimensional Heston
 * equation, on grids with fix Dirichlet-0-Boundaries.
 *
 * @version $HEAD$
 */
class OperationHestonKLinear : public sg::pde::UpDownFourOpDims
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's sg::base::GridStorage object
	 * @param coef vector that contains the constant coefficients of this operation
	 */
	OperationHestonKLinear(sg::base::GridStorage* storage, double***** coef);

	/**
	 * Destructor
	 */
	virtual ~OperationHestonKLinear();

protected:

	// Unidirectional
	void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	// Singles
	void downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	// Doubles
	void downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	// Triples
	void downOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void downOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

	// Quadruples
	void downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);
	void upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim);

};

}
}

#endif /* OPERATIONHESTONKLINEAR_HPP */
