/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/grid/ModWaveletGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg
{
namespace opt
{

ModWaveletGrid::ModWaveletGrid(std::istream &istr) : base::Grid(istr)
{
}

ModWaveletGrid::ModWaveletGrid(size_t dim)
{
    this->storage = new base::GridStorage(dim);
}

ModWaveletGrid::~ModWaveletGrid()
{
}

const char *ModWaveletGrid::getType()
{
    return "modWavelet";
}

base::Grid *ModWaveletGrid::unserialize(std::istream &istr)
{
    return new ModWaveletGrid(istr);
}

base::GridGenerator *ModWaveletGrid::createGridGenerator()
{
    return new base::StandardGridGenerator(this->storage);
}

}
}
