/******************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiDownBBLinearStretched.hpp"

namespace sg
{

namespace detail
{

XPhidPhiDownBBLinearStretched::XPhidPhiDownBBLinearStretched(GridStorage* storage) : storage(storage), stretching(storage->getStretching())
{
}

XPhidPhiDownBBLinearStretched::~XPhidPhiDownBBLinearStretched()
{
}

void XPhidPhiDownBBLinearStretched::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	//	double q = boundingBox->getIntervalWidth(dim);
	//	double t = boundingBox->getIntervalOffset(dim);


	rec(source, result, index, dim, 0.0, 0.0);

}

/*void XPhidPhiDownBBLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{
	size_t seq = index.seq();

	double alpha_value = source[seq];

	GridStorage::index_type::level_type l;
	GridStorage::index_type::index_type i;

	index.get(dim, l, i);

	double hhalf = 1.0/pow(2.0,static_cast<int>(l+1));
	double i_dbl = static_cast<double>(i);

	// integration
	result[seq] = (  ( (fl * ((hhalf * i_dbl) - hhalf)) + (fr * (((-1.0)*(hhalf * i_dbl)) - hhalf)) ) - ((1.0/3.0) * (((1.0/pow(2.0, static_cast<int>(l)))) * alpha_value)) );    // diagonal entry

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

void XPhidPhiDownBBLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
{
	size_t seq = index.seq();

	double alpha_value = source[seq];

	GridStorage::index_type::level_type l;
	GridStorage::index_type::index_type i;

	index.get(dim, l, i);
	//get the positions of the current index as well as its left and right neighbors
	double posl=0, posr=0, posc=0;
//	this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posl, posr );
//	posc = this->stretching->getCoordinates(static_cast<int>(l), static_cast<int>(i), dim);
	this->stretching->getAdjacentPositions(static_cast<int>(l), static_cast<int>(i), dim, posc, posl, posr );
	double baseLength = posr - posl;
	double leftLength = posc - posl;
	double rightLength = posr - posc;

	//	double hhalf = 1.0/pow(2.0,static_cast<int>(l+1));
	//	double i_dbl = static_cast<double>(i);

	// integration
	//	result[seq] = (  ( (fl * ((q*((hhalf * i_dbl) - hhalf))+(0.5*t))) + (fr * ((q*(((-1.0)*(hhalf * i_dbl)) - hhalf))-(0.5*t))) )
	//						  - ((1.0/3.0) * (((1.0/pow(2.0, static_cast<int>(l))) * q) * alpha_value)) );    // diagonal entry

	////Switch the blocks if the result is wrong
	result[seq] = fl*(1.0/6.0)*(2*posc+2*posl-posr) - fr*(1.0/6.0)*(2*posc-posl+2*posr)
								- 1.0/6.0*(baseLength)*alpha_value;	// diagonal entry


	// dehierarchisation
	//	double fm = ((fl+fr)/2.0) + alpha_value;
	double fm =  (fr-fl)*(leftLength)/(baseLength)+fl+ alpha_value;

	if(!index.hint())
	{
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fl, fm);
			//			recBB(source, result, index, dim, fl, fm, q, t);
		}

		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fm, fr);
			//			recBB(source, result, index, dim, fm, fr, q, t);
		}

		index.up(dim);
	}
}

} // namespace detail

} // namespace sg
