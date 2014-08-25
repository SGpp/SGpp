/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/grid/WaveletGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg
{
namespace opt
{

WaveletGrid::WaveletGrid(std::istream &istr) : base::Grid(istr)
{
}

WaveletGrid::WaveletGrid(size_t dim)
{
    this->storage = new base::GridStorage(dim);
}

WaveletGrid::~WaveletGrid()
{
}

const char *WaveletGrid::getType()
{
    return "Wavelet";
}

base::Grid *WaveletGrid::unserialize(std::istream &istr)
{
    return new WaveletGrid(istr);
}

base::GridGenerator *WaveletGrid::createGridGenerator()
{
    return new base::StandardGridGenerator(this->storage);
}

}
}
