/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/BsplineBoundaryGrid.hpp"
#include "base/grid/generation/BoundaryGridGenerator.hpp"
#include "base/exception/factory_exception.hpp"

#include <iostream>

namespace sg
{
namespace base
{

BsplineBoundaryGrid::BsplineBoundaryGrid(std::istream& istr) :
    Grid(istr),
    degree(1 << 16)
{
    istr >> degree;
}


BsplineBoundaryGrid::BsplineBoundaryGrid(size_t dim, size_t degree) : degree(degree)
{
    this->storage = new GridStorage(dim);
}

BsplineBoundaryGrid::~BsplineBoundaryGrid()
{
}

const char *BsplineBoundaryGrid::getType()
{
    return "BsplineBoundary";
}

size_t BsplineBoundaryGrid::getDegree()
{
    return this->degree;
}

Grid *BsplineBoundaryGrid::unserialize(std::istream& istr)
{
    return new BsplineBoundaryGrid(istr);
}

void BsplineBoundaryGrid::serialize(std::ostream& ostr)
{
    this->Grid::serialize(ostr);
    ostr << degree << std::endl;
}

GridGenerator *BsplineBoundaryGrid::createGridGenerator()
{
    return new BoundaryGridGenerator(this->storage);
}

}
}
