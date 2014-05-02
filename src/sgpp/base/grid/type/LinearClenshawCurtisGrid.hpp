/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef LINEARCLENSHAWCURTISGRID_HPP
#define LINEARCLENSHAWCURTISGRID_HPP

#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg
{
namespace base
{

class LinearClenshawCurtisGrid : public Grid
{
protected:
    LinearClenshawCurtisGrid(std::istream& istr);
    
public:
    LinearClenshawCurtisGrid(size_t dim);
    virtual ~LinearClenshawCurtisGrid();
    
    virtual const char* getType();
    
    virtual GridGenerator* createGridGenerator();
    
    static Grid* unserialize(std::istream& istr);
};

}
}

#endif /* LINEARCLENSHAWCURTISGRID_HPP */
