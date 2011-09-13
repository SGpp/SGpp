/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#ifndef OPERATIONEVALMODBSPLINE_HPP
#define OPERATIONEVALMODBSPLINE_HPP

#include "operation/common/OperationEval.hpp"
#include "grid/GridStorage.hpp"
#include "basis/modbspline/modified_bspline_base.hpp"
#include "data/DataVector.hpp"

namespace sg
{
namespace base
{

/**
 * This class implements OperationEval for a grids with modified Bspline basis functions with a certain degree
 *
 * @version $HEAD$
 */
class OperationEvalModBspline : public OperationEval
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 * @param degree the polynom's max. degree
	 */
	OperationEvalModBspline(GridStorage* storage, size_t degree) : storage(storage), base(degree) {}

	/**
	 * Destructor
	 */
	virtual ~OperationEvalModBspline() {}

	virtual double eval(DataVector& alpha, std::vector<double>& point);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
	/// Mod Bspline Basis object
	SModBsplineBase base;
};

}
}

#endif /* OPERATIONEVALMODBSPLINE_HPP */
