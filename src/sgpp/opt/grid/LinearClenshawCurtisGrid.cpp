/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/grid/LinearClenshawCurtisGrid.hpp"
#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include <iostream>

namespace sg
{
namespace opt
{

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(std::istream &istr) : base::Grid(istr)
{
}

LinearClenshawCurtisGrid::LinearClenshawCurtisGrid(size_t dim,
                                                   const tools::CosineTable *cosine_table) :
    cosine_table(cosine_table)
{
    this->storage = new base::GridStorage(dim);
}

LinearClenshawCurtisGrid::~LinearClenshawCurtisGrid()
{
}

const char *LinearClenshawCurtisGrid::getType()
{
    return "linearClenshawCurtis";
}

base::Grid *LinearClenshawCurtisGrid::unserialize(std::istream &istr)
{
    return new LinearClenshawCurtisGrid(istr);
}

base::GridGenerator *LinearClenshawCurtisGrid::createGridGenerator()
{
    return new base::TrapezoidBoundaryGridGenerator(this->storage);
}

const tools::CosineTable *LinearClenshawCurtisGrid::getCosineTable() const
{
    return this->cosine_table;
}

void LinearClenshawCurtisGrid::setCosineTable(const tools::CosineTable *cosine_table)
{
    this->cosine_table = cosine_table;
}

}
}
