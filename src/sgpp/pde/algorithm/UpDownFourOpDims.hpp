/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef UPDOWNFOUROPDIMS_HPP
#define UPDOWNFOUROPDIMS_HPP

#include <vector>

#include "base/grid/GridStorage.hpp"
#include "base/operation/OperationMatrix.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#ifndef TASKS_PARALLEL_UPDOWN
#define TASKS_PARALLEL_UPDOWN 4
#endif

namespace sg
{
namespace pde
{

/**
 * Implements the Up/Down scheme with four dimensions with special operations: i,j,k,l
 *
 * @version $HEAD$
 */
class UpDownFourOpDims: public sg::base::OperationMatrix
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's sg::base::GridStorage object
	 * @param 4d tensor that contains the constant coefficients of this operation
	 */
	UpDownFourOpDims(sg::base::GridStorage* storage, double**** coef);

	/**
	 * Constructor
	 *
	 * @param storage the grid's sg::base::GridStorage object
	 */
	UpDownFourOpDims(sg::base::GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~UpDownFourOpDims();

	virtual void mult(sg::base::DataVector& alpha, sg::base::DataVector& result);

protected:
	typedef sg::base::GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	sg::base::GridStorage* storage;
	/// Pointer to the coefficients of this bilinear form
	double**** coefs;
	/// algorithmic dimensions, operator is applied in this dimensions
	const std::vector<size_t> algoDims;
	/// number of algorithmic dimensions
	const size_t numAlgoDims_;
	/// max number of parallel stages (dimension recursive calls)
	static const size_t maxParallelDims_ = TASKS_PARALLEL_UPDOWN;

	void updown(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);

	void specialOpOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
	void specialOpTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
	void specialOpThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
	void specialOpFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);
	void specialOpX(sg::base::DataVector& alpha, sg::base::DataVector& result, void (sg::pde::UpDownFourOpDims::*pt2UpFunc)(sg::base::DataVector&, sg::base::DataVector&, size_t), void (sg::pde::UpDownFourOpDims::*pt2DownFunc)(sg::base::DataVector&, sg::base::DataVector&, size_t), size_t dim, size_t op_dim_one, size_t op_dim_two, size_t op_dim_three, size_t op_dim_four);


	//	void specialOpOneAndOpTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t op_dim_one, size_t op_dim_two);
	//	virtual void up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	//	virtual void down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	virtual void downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	virtual void upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	virtual void downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	virtual void upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	virtual void downOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	virtual void upOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	virtual void downOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	virtual void upOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	//	virtual void downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
	//	virtual void upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) = 0;
};

}
}

#endif /* UPDOWNFOUROPDIMS_HPP */
