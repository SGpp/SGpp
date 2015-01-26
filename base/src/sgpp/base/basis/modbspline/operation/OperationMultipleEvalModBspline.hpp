/* ****************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)
#ifndef OPERATIONMULTIPLEEVALMODBSPLINE_HPP
#define OPERATIONMULTIPLEEVALMODBSPLINE_HPP

#include <sgpp/base/operation/OperationMultipleEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/basis/modbspline/ModifiedBsplineBasis.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class implements OperationMultipleEval for a grid with modified Bspline basis functions
 *
 * @version $HEAD$
 */
class OperationMultipleEvalModBspline: public OperationMultipleEval {
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the Bspline's degree
	 * @param dataset the dataset that should be evaluated
	 */
	OperationMultipleEvalModBspline(Grid &grid, size_t degree,
			DataMatrix &dataset) :
			OperationMultipleEval(grid, dataset), base(degree) {
		this->storage = grid.getStorage();
	}

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalModBspline() {
	}

	virtual void mult(DataVector& alpha, DataVector& result);
	virtual void multTranspose(DataVector& source, DataVector& result);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Mod Bspline Basis object
	SModBsplineBase base;
};

}
}

#endif /* OPERATIONMULTIPLEEVALMODBSPLINE_HPP */
