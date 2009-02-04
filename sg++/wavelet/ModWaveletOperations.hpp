/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2008 Deepak Pandey

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MODWAVELETOPERATIONS_HPP_
#define MODWAVELETOPERATIONS_HPP_

#include "Operations.hpp"

namespace sg
{

class OperationBModWavelet : public OperationB
{
public:
	OperationBModWavelet(GridStorage* storage) : storage(storage) {}
	virtual ~OperationBModWavelet() {}

	virtual void mult(DataVector& alpha, DataVector& data, DataVector& result);
	virtual void multTranspose(DataVector& alpha, DataVector& data, DataVector& result);

protected:
	GridStorage* storage;
};


class OperationEvalModWavelet : public OperationEval
{
public:
	OperationEvalModWavelet(GridStorage* storage) : storage(storage) {}
	virtual ~OperationEvalModWavelet() {}

	virtual double eval(DataVector& alpha, std::vector<double>& point);
	virtual double test(DataVector& alpha, DataVector& data, DataVector& classes);

protected:
	GridStorage* storage;

};

/**
 * Hierarchisation on sparse grid, mod wavelet case
 */
class OperationHierarchisationModWavelet : public OperationHierarchisation
{
public:
	OperationHierarchisationModWavelet(GridStorage* storage) : storage(storage) {}
	virtual ~OperationHierarchisationModWavelet() {}

	virtual void doHierarchisation(DataVector& node_values);
	virtual void doDehierarchisation(DataVector& alpha);

protected:
	GridStorage* storage;
};

}

#endif /*MODWAVELETOPERATIONS_HPP_*/
