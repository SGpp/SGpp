/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef BSPLINECLENSHAWCURTISGRID_HPP
#define BSPLINECLENSHAWCURTISGRID_HPP

#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg
{
namespace base
{

class BsplineClenshawCurtisGrid : public Grid
{
protected:
    BsplineClenshawCurtisGrid(std::istream& istr);
    
public:
    BsplineClenshawCurtisGrid(size_t dim, size_t degree);
    virtual ~BsplineClenshawCurtisGrid();
    
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

#endif /* BSPLINECLENSHAWCURTISGRID_HPP */
