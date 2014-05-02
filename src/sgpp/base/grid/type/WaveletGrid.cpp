/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/WaveletGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"
#include "base/exception/factory_exception.hpp"

#include <iostream>

namespace sg
{
namespace base
{

WaveletGrid::WaveletGrid(std::istream& istr) : Grid(istr)
{
}

WaveletGrid::WaveletGrid(size_t dim)
{
    this->storage = new GridStorage(dim);
}

WaveletGrid::~WaveletGrid()
{
}

const char *WaveletGrid::getType()
{
    return "Wavelet";
}

Grid *WaveletGrid::unserialize(std::istream& istr)
{
    return new WaveletGrid(istr);
}

GridGenerator *WaveletGrid::createGridGenerator()
{
    return new StandardGridGenerator(this->storage);
}

}
}
