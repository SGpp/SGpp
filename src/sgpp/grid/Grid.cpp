/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/Grid.hpp"
#include "grid/type/LinearGrid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"
#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "grid/type/ModLinearGrid.hpp"

#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "operation/common/OperationIdentity.hpp"

#include "exception/factory_exception.hpp"

#include <iostream>
#include <sstream>

namespace sg
{

Grid* Grid::createLinearGrid(size_t dim)
{
	return new LinearGrid(dim);
}

Grid* Grid::createLinearBoundaryGrid(size_t dim)
{
	return new LinearBoundaryGrid(dim);
}

Grid* Grid::createLinearTrapezoidBoundaryGrid(size_t dim)
{
	return new LinearTrapezoidBoundaryGrid(dim);
}

Grid* Grid::createModLinearGrid(size_t dim)
{
	return new ModLinearGrid(dim);
}

Grid* Grid::createPolyGrid(size_t dim, size_t degree)
{
  //	return new PolyGrid(dim, degree);
	throw factory_exception("Grid-type not in this revision");
	return NULL;
}

Grid* Grid::createModWaveletGrid(size_t dim)
{
  //    return new ModWaveletGrid(dim);
	throw factory_exception("Grid-type not in this revision");
	return NULL;
}

Grid* Grid::createModBsplineGrid(size_t dim, size_t degree)
{
  //    return new ModBsplineGrid(dim, degree);
	throw factory_exception("Grid-type not in this revision");
	return NULL;
}

OperationMatrix* Grid::createOperationIdentity()
{
	return new OperationIdentity();
}

Grid* Grid::createModPolyGrid(size_t dim, size_t degree)
{
  //	return new ModPolyGrid(dim, degree);
	throw factory_exception("Grid-type not in this revision");
	return NULL;
}

Grid* Grid::unserialize(const std::string& istr)
{
	std::istringstream istream;
	istream.str(istr);

	return Grid::unserialize(istream);
}

Grid* Grid::unserialize(std::istream& istr)
{
	std::string gridtype;
	istr >> gridtype;

	if(typeMap().count(gridtype) > 0)
	{
		// This is pretty esoteric. It calls a function pointer out of a map.
		return typeMap()[gridtype](istr);
	}
	else
	{
		throw factory_exception("factory_exception unserialize: unknown gridtype");
	}

	return NULL;
}

std::map<std::string, Grid::Factory>& Grid::typeMap()
{
	// This is only executed once!
	static factoryMap* tMap = new factoryMap();
	if(tMap->size() == 0)
	{
		/*
		 * Insert factories here. This methods may NOT read the grid type.
		 * This map takes a string, function pointer pair.
		 */
		tMap->insert(std::make_pair("NULL",Grid::nullFactory));
		tMap->insert(std::make_pair("linear", LinearGrid::unserialize));
		tMap->insert(std::make_pair("linearBoundary", LinearBoundaryGrid::unserialize));
		tMap->insert(std::make_pair("linearTrapezoidBoundary", LinearTrapezoidBoundaryGrid::unserialize));
		tMap->insert(std::make_pair("modlinear", ModLinearGrid::unserialize));
		/*
		tMap->insert(std::make_pair("poly", PolyGrid::unserialize));
		tMap->insert(std::make_pair("modpoly", ModPolyGrid::unserialize));
		tMap->insert(std::make_pair("modWavelet", ModWaveletGrid::unserialize));
		tMap->insert(std::make_pair("modBspline", ModBsplineGrid::unserialize));
		*/
	}

	return *tMap;
}

/**
 * Factory for everything we don't know.
 */
Grid* Grid::nullFactory(std::istream&)
{
	throw factory_exception("factory_exeception unserialize: unsupported gridtype");
	return NULL;
}

Grid::Grid(std::istream& istr) : storage(NULL)
{
	int hasStorage;
	istr >> hasStorage;
	if(hasStorage == 1)
	{
		storage = new GridStorage(istr);
	}
}

Grid::Grid() : storage(NULL)
{
}

Grid::~Grid()
{
	if(storage != NULL)
	{
		delete storage;
	}

	if(evalOp != NULL)
	{
		delete evalOp;
		evalOp = NULL;
	}
}

GridStorage* Grid::getStorage()
{
	return this->storage;
}

BoundingBox* Grid::getBoundingBox()
{
	return this->storage->getBoundingBox();
}

void Grid::setBoundingBox(BoundingBox& bb)
{
	this->storage->setBoundingBox(bb);
}

void Grid::serialize(std::string& ostr)
{
	std::ostringstream ostream;
	this->serialize(ostream);

	ostr = ostream.str();
}

std::string Grid::serialize()
{
	std::ostringstream ostream;
	this->serialize(ostream);

	return ostream.str();
}

void Grid::serialize(std::ostream& ostr)
{
	ostr << this->getType() << std::endl;
	if(storage != NULL)
	{
		ostr << "1" << std::endl;
		storage->serialize(ostr);
	}
	else
	{
		ostr << "0" << std::endl;
	}
}

void Grid::refine(DataVector* vector, int numOfPoints)
{
	// @todo (khakhutv) (low) different refinemente Functors
	this->createGridGenerator()->refine(new SurplusRefinementFunctor(vector, numOfPoints));
}

OperationEval* Grid::evalOp(NULL);

double Grid::eval(DataVector& alpha, DataVector& point){
	if(this->evalOp == NULL) this->evalOp = this->createOperationEval();
	return this->evalOp->eval(alpha, point);
}

void Grid::insertPoint(size_t dim, unsigned int levels[], unsigned int indices[], bool isLeaf){
	//create HashGridIndex object for the point
	GridIndex pointIndex = new GridIndex(dim);
	for (unsigned int i=0; i<dim-1; i++){
		pointIndex.push(i, levels[i], indices[i]);
	}
	//insert last level/index and hash
	pointIndex.set(dim-1, levels[dim-1], indices[dim-1], isLeaf);
	//insert point to the GridStorage
	storage->insert(pointIndex);
}

int Grid::getSize(){
	return this->storage->size();
}

}
