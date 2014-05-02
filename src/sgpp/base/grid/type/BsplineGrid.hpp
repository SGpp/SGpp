/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef BSPLINEGRID_HPP
#define BSPLINEGRID_HPP

#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg
{
namespace base
{

class BsplineGrid : public Grid
{
protected:
    BsplineGrid(std::istream& istr);
    
public:
    BsplineGrid(size_t dim, size_t degree);
    virtual ~BsplineGrid();
    
    virtual const char *getType();
    
    virtual GridGenerator *createGridGenerator();
    
    static Grid *unserialize(std::istream& istr);
    
    virtual void serialize(std::ostream& ostr);
    virtual size_t getDegree();
    
protected:
    size_t degree;
};

}
}

#endif /* BSPLINEGRID_HPP */
