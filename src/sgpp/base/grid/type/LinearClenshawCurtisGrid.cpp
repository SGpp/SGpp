/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/LinearClenshawCurtisGrid.hpp"
#include "base/grid/generation/BoundaryGridGenerator.hpp"
#include "base/exception/factory_exception.hpp"

#include <iostream>

namespace sg
{
namespace base
{

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(std::istream& istr) : Grid(istr)
{
}

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(size_t dim)
{
    this->storage = new GridStorage(dim);
}

LinearClenshawCurtisGrid::~LinearClenshawCurtisGrid()
{
}

const char* LinearClenshawCurtisGrid::getType()
{
    return "linearClenshawCurtis";
}

Grid* LinearClenshawCurtisGrid::unserialize(std::istream& istr)
{
    return new LinearClenshawCurtisGrid(istr);
}

GridGenerator* LinearClenshawCurtisGrid::createGridGenerator()
{
    return new BoundaryGridGenerator(this->storage);
}

}
}
