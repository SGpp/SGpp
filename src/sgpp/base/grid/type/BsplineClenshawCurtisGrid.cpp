/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/BsplineClenshawCurtisGrid.hpp"
#include "base/grid/generation/BoundaryGridGenerator.hpp"
#include "base/exception/factory_exception.hpp"

#include <iostream>

namespace sg
{
namespace base
{

// sg::base::GridStorageClenshawCurtis* {aka
// sg::base::HashGridStorage<sg::base::HashGridIndexClenshawCurtis<unsigned int, unsigned int> >*}
// sg::base::GridStorage* {aka
// sg::base::HashGridStorage<sg::base::HashGridIndex<unsigned int, unsigned int> >*}

/*class GridStorageClenshawCurtis : public GridStorage
{
}*/
// TODO (can be removed)

BsplineClenshawCurtisGrid::BsplineClenshawCurtisGrid(std::istream& istr) :
    Grid(istr),
    degree(1 << 16)
{
    istr >> degree;
}

//template class HashGridStorage<GridIndexClenshawCurtis>;

BsplineClenshawCurtisGrid::BsplineClenshawCurtisGrid(size_t dim, size_t degree) : degree(degree)
{
    this->storage = new GridStorage(dim);
}

BsplineClenshawCurtisGrid::~BsplineClenshawCurtisGrid()
{
}

const char *BsplineClenshawCurtisGrid::getType()
{
    return "BsplineClenshawCurtis";
}

size_t BsplineClenshawCurtisGrid::getDegree()
{
    return this->degree;
}

Grid *BsplineClenshawCurtisGrid::unserialize(std::istream& istr)
{
    return new BsplineClenshawCurtisGrid(istr);
}

void BsplineClenshawCurtisGrid::serialize(std::ostream& ostr)
{
    this->Grid::serialize(ostr);
    ostr << degree << std::endl;
}

GridGenerator *BsplineClenshawCurtisGrid::createGridGenerator()
{
    return new BoundaryGridGenerator(this->storage);
}

}
}
