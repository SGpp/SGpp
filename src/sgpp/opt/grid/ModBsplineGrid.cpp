/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/grid/ModBsplineGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg
{
namespace opt
{

ModBsplineGrid::ModBsplineGrid(std::istream &istr) : base::Grid(istr), degree(1 << 16)
{
    istr >> degree;
}

ModBsplineGrid::ModBsplineGrid(size_t dim, size_t degree) : degree(degree)
{
    this->storage = new base::GridStorage(dim);
}

ModBsplineGrid::~ModBsplineGrid()
{
}

const char *ModBsplineGrid::getType()
{
    return "modBspline";
}

size_t ModBsplineGrid::getDegree()
{
    return this->degree;
}

base::Grid *ModBsplineGrid::unserialize(std::istream &istr)
{
    return new ModBsplineGrid(istr);
}

void ModBsplineGrid::serialize(std::ostream &ostr)
{
    this->base::Grid::serialize(ostr);
    ostr << degree << std::endl;
}

base::GridGenerator *ModBsplineGrid::createGridGenerator()
{
    return new base::StandardGridGenerator(this->storage);
}

}
}
