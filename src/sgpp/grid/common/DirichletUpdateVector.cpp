/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/common/DirichletUpdateVector.hpp"

namespace sg
{

DirichletUpdateVector::DirichletUpdateVector(GridStorage* storage): myBoundingBox(storage->getBoundingBox()), storage(storage)
{
}

DirichletUpdateVector::~DirichletUpdateVector()
{
}

void DirichletUpdateVector::applyDirichletConditions(DataVector& updateVector, DataVector& sourceVector)
{
	for (size_t i = 0; i < storage->size(); i++)
	{
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == false)
		{
			updateVector.set(i, sourceVector.get(i));
		}
	}
}

void DirichletUpdateVector::setBoundariesToZero(DataVector& updateVector)
{
	for (size_t i = 0; i < storage->size(); i++)
	{
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == false)
		{
			updateVector.set(i, 0.0);
		}
	}
}

void DirichletUpdateVector::setInnerPointsToZero(DataVector& updateVector)
{
	for (size_t i = 0; i < storage->size(); i++)
	{
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == true)
		{
			updateVector.set(i, 0.0);
		}
	}
}

void DirichletUpdateVector::multiplyBoundary(DataVector& updateVector, double value)
{
	for (size_t i = 0; i < storage->size(); i++)
	{
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == false)
		{
			updateVector.set(i, updateVector.get(i)*value);
		}
	}
}

void DirichletUpdateVector::getfactor(DataVector& factor, double T)
{
	double tmp;
	for (size_t i = 0; i < storage->size(); i++)
	{
		std::string coords = (*storage)[i]->getCoordsStringBB(*this->myBoundingBox);
		std::stringstream coordsStream(coords);
		double* dblFuncValues = new double[2];
		for (size_t j = 0; j < 2; j++)
		{
			coordsStream >> tmp;
			dblFuncValues[j] = tmp;
		}
		//std::cout<<dblFuncValues[1]<<std::endl;
		factor.set(i, exp((-1.0)*dblFuncValues[1]*T));
	}
}

void DirichletUpdateVector::multiplyBoundaryVector(DataVector& updateVector,DataVector& factor)
{
	for (size_t i = 0; i < storage->size(); i++)
	{
		GridIndex* curPoint = (*storage)[i];
		if (curPoint->isInnerPoint() == false)
		{
			updateVector.set(i, updateVector.get(i)* factor.get(i));
		}
	}
}

void DirichletUpdateVector::multiplyrBSHW(DataVector& updateVector)
{
	double tmp;
	for (size_t i = 0; i < storage->size(); i++)
	{
		std::string coords = (*storage)[i]->getCoordsStringBB(*this->myBoundingBox);
		std::stringstream coordsStream(coords);
		double* dblFuncValues = new double[2];
        for (size_t j = 0; j < 2; j++)
			{
			    coordsStream >> tmp;
                dblFuncValues[j] = tmp;
			}
       // std::cout<< dblFuncValues[1]<< std::endl;
		updateVector.set(i, updateVector.get(i)* dblFuncValues[1]);
	}
}
}
