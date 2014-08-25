/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/grid/WaveletTrapezoidBoundaryGrid.hpp"
#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"

namespace sg
{
namespace opt
{

WaveletTrapezoidBoundaryGrid::WaveletTrapezoidBoundaryGrid(std::istream &istr) :
    base::Grid(istr)
{
}

WaveletTrapezoidBoundaryGrid::WaveletTrapezoidBoundaryGrid(size_t dim)
{
    this->storage = new base::GridStorage(dim);
}

WaveletTrapezoidBoundaryGrid::~WaveletTrapezoidBoundaryGrid()
{
}

const char *WaveletTrapezoidBoundaryGrid::getType()
{
    return "WaveletTrapezoidBoundary";
}

base::Grid *WaveletTrapezoidBoundaryGrid::unserialize(std::istream &istr)
{
    return new WaveletTrapezoidBoundaryGrid(istr);
}

base::GridGenerator *WaveletTrapezoidBoundaryGrid::createGridGenerator()
{
    return new base::TrapezoidBoundaryGridGenerator(this->storage);
}

}
}
