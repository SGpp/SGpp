/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/WaveletBoundaryGrid.hpp"
#include "base/grid/generation/BoundaryGridGenerator.hpp"
#include "base/exception/factory_exception.hpp"

#include <iostream>

namespace sg
{
namespace base
{

WaveletBoundaryGrid::WaveletBoundaryGrid(std::istream& istr) : Grid(istr)
{
}

WaveletBoundaryGrid::WaveletBoundaryGrid(size_t dim)
{
    this->storage = new GridStorage(dim);
}

WaveletBoundaryGrid::~WaveletBoundaryGrid()
{
}

const char *WaveletBoundaryGrid::getType()
{
    return "WaveletBoundary";
}

Grid *WaveletBoundaryGrid::unserialize(std::istream& istr)
{
    return new WaveletBoundaryGrid(istr);
}

GridGenerator *WaveletBoundaryGrid::createGridGenerator()
{
    return new BoundaryGridGenerator(this->storage);
}

}
}
