// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALMODPOLY_HPP
#define OPERATIONMULTIPLEEVALMODPOLY_HPP

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationMultipleEval for a grids with mod poly basis ansatzfunctions
 *
 * @version $HEAD$
 */
class OperationMultipleEvalModPoly: public OperationMultipleEval {
public:
	/**
	 * Constructor
	 *
	 * @param grid grid
	 * @param degree the polynom's max. degree
	 * @param dataset the dataset that should be evaluated
	 */
	OperationMultipleEvalModPoly(Grid &grid, size_t degree, DataMatrix &dataset) :
			OperationMultipleEval(grid, dataset), base(degree) {
		this->storage = grid.getStorage();
	}

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalModPoly() {
	}

	virtual void mult(DataVector& alpha, DataVector& result);
	virtual void multTranspose(DataVector& source, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Mod Poly Basis object
	SPolyModifiedBase base;
};

}
}

#endif /* OPERATIONMULTIPLEEVALMODPOLY_HPP */
