/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp"

namespace sg
{

namespace detail
{

PhiPhiDownBBLinearStretched::PhiPhiDownBBLinearStretched(GridStorage* storage) : storage(storage), stretching(storage->getStretching())
{
}

PhiPhiDownBBLinearStretched::~PhiPhiDownBBLinearStretched()
{
}

void PhiPhiDownBBLinearStretched::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
		rec(source, result, index, dim, 0.0, 0.0);

}

void PhiPhiDownBBLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{

	size_t seq = index.seq();

	double alpha_value = source[seq];

	GridStorage::index_type::level_type current_level;
	GridStorage::index_type::index_type current_index;

	index.get(dim, current_level, current_index);

	double posl=0, posr=0, currentPosition=0;
//	this->stretching->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, posl, posr );
//	currentPosition = this->stretching->getCoordinates(static_cast<int>(current_level), static_cast<int>(current_index), dim);
	this->stretching->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, currentPosition,posl, posr );
	double baseLength = posr - posl;
	double leftLength = currentPosition - posl;
	double rightLength = posr - currentPosition;
//	double h = (1.0/(pow(2.0, static_cast<int>(l))));

	// integration
	result[seq] = (1.0/3.0)*(baseLength)*alpha_value + fl/6.0*(baseLength+rightLength) + fr/6.0*(baseLength+leftLength);
//	result[seq] = (h * ((fl+fr)/2.0)) + (2.0/3.0)*h * alpha_value);

	// dehierarchisation
//	double fm = ((fl+fr)/2.0) + alpha_value;
	double fm  = (fr-fl)*(leftLength)/(baseLength)+fl + alpha_value;


	if(!index.hint())
	{
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fl, fm);
		}

		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fm, fr);
		}

		index.up(dim);
	}
}

/*void PhiPhiDownBBLinearStretched::recBB(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr, double q, double t)
{
	size_t seq = index.seq();

	double alpha_value = source[seq];

	GridStorage::index_type::level_type l;
	GridStorage::index_type::index_type i;

	index.get(dim, l, i);

	double h = (1.0/(pow(2.0, static_cast<int>(l))));

	// integration
	result[seq] = ((h * ((fl+fr)/2.0)) * q) + ((((2.0/3.0) * h) * alpha_value) * q);    // diagonal entry

	// dehierarchisation
	double fm = ((fl+fr)/2.0) + alpha_value;

	if(!index.hint())
	{
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			recBB(source, result, index, dim, fl, fm, q, t);
		}

		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			recBB(source, result, index, dim, fm, fr, q, t);
		}

		index.up(dim);
	}
}*/

} // namespace detail

} // namespace sg
