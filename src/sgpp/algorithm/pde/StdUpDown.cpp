/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "algorithm/pde/StdUpDown.hpp"
using namespace sg::base;

namespace sg
{

StdUpDown::StdUpDown(GridStorage* storage)
{
	this->storage = storage;
	this->algoDims = this->storage->getAlgorithmicDimensions();
}

StdUpDown::~StdUpDown()
{
}

void StdUpDown::mult(DataVector& alpha, DataVector& result)
{
	DataVector beta(result.getSize());
	result.setAll(0.0);
	#pragma omp parallel
	{
		#pragma omp single nowait
		{
			this->updown(alpha, beta, this->algoDims.size() - 1);
		}
	}
	result.add(beta);
}

void StdUpDown::multParallelBuildingBlock(DataVector& alpha, DataVector& result)
{
	DataVector beta(result.getSize());
	result.setAll(0.0);

	this->updown(alpha, beta, this->algoDims.size() - 1);

	result.add(beta);
}

void StdUpDown::updown(DataVector& alpha, DataVector& result, size_t dim)
{
	//Unidirectional scheme
	if(dim > 0)
	{
		// Reordering ups and downs
		DataVector temp(alpha.getSize());
		DataVector result_temp(alpha.getSize());
		DataVector temp_two(alpha.getSize());

		#pragma omp task shared(alpha, temp, result)
		{
			up(alpha, temp, this->algoDims[dim]);
			updown(temp, result, dim-1);
		}


		#pragma omp task shared(alpha, temp_two, result_temp)
		{
			updown(alpha, temp_two, dim-1);
			down(temp_two, result_temp, this->algoDims[dim]);
		}

		#pragma omp taskwait

		result.add(result_temp);
	}
	else
	{
		// Terminates dimension recursion
		DataVector temp(alpha.getSize());

		#pragma omp task shared(alpha, result)
		up(alpha, result, this->algoDims[dim]);

		#pragma omp task shared(alpha, temp)
		down(alpha, temp, this->algoDims[dim]);

		#pragma omp taskwait

		result.add(temp);
	}
}

}
