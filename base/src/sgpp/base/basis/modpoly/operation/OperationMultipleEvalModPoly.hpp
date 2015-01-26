/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)
#ifndef OPERATIONMULTIPLEEVALMODPOLY_HPP
#define OPERATIONMULTIPLEEVALMODPOLY_HPP

#include <sgpp/base/operation/OperationMultipleEval.hpp>
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
	 * @param storage the grid's GridStorage object
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
	SModPolyBase base;
};

}
}

#endif /* OPERATIONMULTIPLEEVALMODPOLY_HPP */
