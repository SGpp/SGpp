/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de, Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "basis/linearstretched/noboundary/common/DowndPhidPhiBBIterativeLinearStretched.hpp"
#include "grid/common/Stretching.hpp"
using namespace sg::base;

namespace sg
{
namespace pde
{

DowndPhidPhiBBIterativeLinearStretched::DowndPhidPhiBBIterativeLinearStretched(GridStorage* storage) : storage(storage)
{
}

DowndPhidPhiBBIterativeLinearStretched::~DowndPhidPhiBBIterativeLinearStretched()
{
}

void DowndPhidPhiBBIterativeLinearStretched::operator()(DataVector& alpha, DataVector& result, size_t dim)
{
	Stretching* stretching= this->storage->getStretching();
//	double q = stretching->getIntervalWidth(dim);
//
//	double Qqout = 1.0/q;

/*	// init the coefficients of the ansatz functions with boundary
	result.setAll(0.0);

	if (q != 1.0)
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			//only affects the diagonal of the stiffness matrix
			result[i] = alpha[i]*(Qqout*pow(2.0, static_cast<int>(level+1)));
		}
	}
	else
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			//only affects the diagonal of the stiffness matrix
			result[i] = alpha[i]*pow(2.0, static_cast<int>(level+1));
		}
	}*/
//	else{
//		std::cout<<"else called"<<std::endl;
//	}
	result.setAll(0.0);

	// traverse all basis function by sequence number
			for(size_t i = 0; i < storage->size(); i++)
			{
				GridStorage::index_type::level_type level;
				GridStorage::index_type::index_type index;
				(*storage)[i]->get(dim, level, index);
				double posl=0, posr=0, posc=0;
				stretching->getAdjacentPositions(static_cast<int>(level), static_cast<int>(index), dim, posc, posl, posr );
				double baseLength = posr - posl;
				double leftLength = posc - posl;
				double rightLength = posr - posc;
				//only affects the diagonal of the stiffness matrix
				result[i] = alpha[i]*baseLength/(leftLength*rightLength);
			}

}

}
}
