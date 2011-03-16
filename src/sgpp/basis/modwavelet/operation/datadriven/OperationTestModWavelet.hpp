/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (dirk.pflueger@in.tum.de)

#ifndef OPERATIONTESTMODWAVELET_HPP
#define OPERATIONTESTMODWAVELET_HPP

#include "operation/datadriven/OperationTest.hpp"
#include "grid/GridStorage.hpp"
using namespace sg::base;

namespace sg
{

/**
 * This class implements OperationTest for a grid with mod wavelet basis ansatzfunctions
 *
 * @version $HEAD$
 */
class OperationTestModWavelet : public OperationTest
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's GridStorage object
	 */
	OperationTestModWavelet(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~OperationTestModWavelet() {}

	virtual double test(DataVector& alpha, DataMatrix& data, DataVector& classes);
	virtual double testMSE(DataVector& alpha, DataMatrix& data, DataVector& refValues);
	virtual double testWithCharacteristicNumber(DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers);

protected:
	/// Pointer to GridStorage object
	GridStorage* storage;
};

}

#endif /* OPERATIONTESTMODWAVELET_HPP */
