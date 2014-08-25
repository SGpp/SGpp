/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/grid/ModLinearGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg
{
namespace opt
{

ModLinearGrid::ModLinearGrid(std::istream &istr) : base::Grid(istr)
{
}

ModLinearGrid::ModLinearGrid(size_t dim)
{
    this->storage = new base::GridStorage(dim);
}

ModLinearGrid::~ModLinearGrid()
{
}

const char *ModLinearGrid::getType()
{
    return "modLinear";
}

base::Grid *ModLinearGrid::unserialize(std::istream &istr)
{
    return new ModLinearGrid(istr);
}

base::GridGenerator *ModLinearGrid::createGridGenerator()
{
    return new base::StandardGridGenerator(this->storage);
}

}
}
