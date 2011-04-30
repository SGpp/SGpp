/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#ifndef OPERATIONTESTLINEARSTRETCHED_HPP
#define OPERATIONTESTLINEARSTRETCHED_HPP

#include "operation/datadriven/OperationTest.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{
namespace datadriven
{

/**
 * This class implements OperationTest for a grids with linearstretched basis ansatzfunctions without boundaries
 */
class OperationTestLinearStretched : public OperationTest
{
public:
	/**
	 * Construtor of OperationTestLinearStretched
	 *
	 * @param storage Pointer to the grid's gridstorage obejct
	 */
	OperationTestLinearStretched(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestLinearStretched() {}

	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes);
	virtual double testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues);
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers);

protected:
	/// Pointer to the grid's gridstorage object
	GridStorage* storage;
};

}
}

#endif /* OPERATIONTESTLINEARSTRETCHED_HPP */
