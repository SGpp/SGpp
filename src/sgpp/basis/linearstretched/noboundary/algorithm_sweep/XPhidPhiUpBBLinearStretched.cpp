/******************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiUpBBLinearStretched.hpp"

namespace sg
{

namespace detail
{

XPhidPhiUpBBLinearStretched::XPhidPhiUpBBLinearStretched(GridStorage* storage) : storage(storage), stretching(storage->getStretching())
{
}

XPhidPhiUpBBLinearStretched::~XPhidPhiUpBBLinearStretched()
{
}

void XPhidPhiUpBBLinearStretched::operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
{
	//	double q = boundingBox->getIntervalWidth(dim);
	//	double t = boundingBox->getIntervalOffset(dim);



	// get boundary values
	double fl = 0.0;
	double fr = 0.0;

	rec(source, result, index, dim, fl, fr);

}

/*void XPhidPhiUpBBLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr)
{
	size_t seq = index.seq();

	fl = fr = 0.0;
	double fml = 0.0;
	double fmr = 0.0;

	GridStorage::index_type::level_type current_level;
	GridStorage::index_type::index_type current_index;

	if(!index.hint())
	{
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fl, fml);
		}

		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fmr, fr);
		}

		index.up(dim);
	}

	index.get(dim, current_level, current_index);

	double fm = fml + fmr;

	double alpha_value = source[seq];

	double helper = (1.0/pow(2.0, static_cast<int>(current_level+1))) * (static_cast<double>(current_index));

	// transposed operations:
	result[seq] = fm;

	fl = (fm/2.0) + ((alpha_value*((-1.0)*helper)) + fl);
	fr = (fm/2.0) + ((alpha_value*helper) + fr);
}*/

void XPhidPhiUpBBLinearStretched::rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr)
{
	size_t seq = index.seq();

	fl = fr = 0.0;
	double fml = 0.0;
	double fmr = 0.0;

	GridStorage::index_type::level_type current_level;
	GridStorage::index_type::index_type current_index;

	if(!index.hint())
	{
		index.left_child(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fl, fml);
			//			recBB(source, result, index, dim, fl, fml, q, t);
		}

		index.step_right(dim);
		if(!storage->end(index.seq()))
		{
			rec(source, result, index, dim, fmr, fr);
			//			recBB(source, result, index, dim, fmr, fr, q, t);
		}

		index.up(dim);
	}

	index.get(dim, current_level, current_index);
	//get the positions of the current index as well as its left and right neighbors
	double posl=0, posr=0, posc=0;
//	this->stretching->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, posl, posr );
//	posc = this->stretching->getCoordinates(static_cast<int>(current_level), static_cast<int>(current_index), dim);
	this->stretching->getAdjacentPositions(static_cast<int>(current_level), static_cast<int>(current_index), dim, posc, posl, posr );
	double baseLength = posr - posl;
	double leftLength = posc - posl;
	double rightLength = posr - posc;

	double fm = fml + fmr;

	double alpha_value = source[seq];

	//	double helper = (1.0/pow(2.0, static_cast<int>(current_level+1))) * (q * static_cast<double>(current_index)) + (0.5*t);

	// transposed operations:
	result[seq] = fm;
	//	fl = (fm/2.0) + ((alpha_value*((-1.0)*helper)) + fl);
	//	fr = (fm/2.0) + ((alpha_value*helper) + fr);

	////Switch the blocks if the result is wrong
//	fl = (1.0/6.0)*(-1*alpha_value*(baseLength+leftLength)) + fl + fm*(leftLength/baseLength);
//	fr = (1.0/6.0)*(-1*alpha_value*(baseLength+rightLength)) + fr + fm*(rightLength/baseLength);
	double c = -1.0/6.0*(posc+posl+posr);
	fl = (c*alpha_value) + fl + fm*(rightLength/baseLength);
	fr = (-c*alpha_value) + fr + fm*(leftLength/baseLength);
//	fl = ((1.0/6.0)*(2*posc+2*posl-posr)*alpha_value) + fl + fm*(leftLength/baseLength);
//	fr = (-(1.0/6.0)*(2*posc-posl+2*posr)*alpha_value) + fr + fm*(rightLength/baseLength);
}

} // namespace detail

} // namespace sg
