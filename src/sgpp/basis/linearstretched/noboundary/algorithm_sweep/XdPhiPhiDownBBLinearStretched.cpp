/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiDownBBLinearStretched.hpp"

namespace sg
{

namespace detail
{

XdPhiPhiDownBBLinearStretched::XdPhiPhiDownBBLinearStretched(GridStorage* storage) : storage(storage), stretching(storage->getStretching())
{
}

XdPhiPhiDownBBLinearStretched::~XdPhiPhiDownBBLinearStretched()
{
}

void XdPhiPhiDownBBLinearStretched::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
//	double q = boundingBox->getIntervalWidth(dim);
//	double t = boundingBox->getIntervalOffset(dim);
//
//
//
//	if (useBB)
//	{
//		recBB(source, result, index, dim, 0.0, 0.0, q, t);
//	}
//	else
//	{
		rec(source, result, index, dim, 0.0, 0.0);
//	}
}

/*void XdPhiPhiDownBBLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{
	//TODO: Fix this (selcuk)
	size_t seq = index.seq();

	double alpha_value = source[seq];

	GridStorage::index_type::level_type l;
	GridStorage::index_type::index_type i;

	index.get(dim, l, i);

	double helper = (1.0/pow(2.0, static_cast<int>(l+1))) * (static_cast<double>(i));

	// integration
	result[seq] = (  ( (fr-fl) * (helper) ) - ((1.0/3.0) * (((1.0/pow(2.0, static_cast<int>(l)))) * alpha_value)) );    // diagonal entry

	// dehierarchisation
	double fm = ((fl+fr)/2.0) + alpha_value;

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
}*/

void XdPhiPhiDownBBLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{
	size_t seq = index.seq();

	double alpha_value = source[seq];

	GridStorage::index_type::level_type l;
	GridStorage::index_type::index_type i;

	index.get(dim, l, i);
	//get the positions of the current index as well as its left and right neighbors
	double posl=0, posr=0, posc=0;
	this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posc, posl, posr );

	double baseLength = posr - posl;
	double leftLength = posc - posl;
	double rightLength = posr - posc;

//	double helper = (1.0/pow(2.0, static_cast<int>(l+1))) * (q * static_cast<double>(i));

	// integration
//	result[seq] = (  ( (fr-fl) * (helper + (0.5*t)) )
//						  - ((1.0/3.0) * (((1.0/pow(2.0, static_cast<int>(l))) * q) * alpha_value)) );    // diagonal entry

	result[seq] = fl*(-1.0/6.0)*(posc+posl+posr)-fr*(-1.0/6.0)*(posc+posl+posr)
					- (1.0/6.0)*baseLength*alpha_value;	// diagonal entry




	// dehierarchisation
//	double fm = ((fl+fr)/2.0) + alpha_value;
	double fm =  (fr-fl)*(leftLength)/(baseLength)+fl+ alpha_value;


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

} // namespace detail

} // namespace sg
