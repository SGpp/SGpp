/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef WAVELETGRID_HPP
#define WAVELETGRID_HPP

#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg
{
namespace base
{

class WaveletGrid : public Grid
{
protected:
    WaveletGrid(std::istream& istr);
    
public:
    WaveletGrid(size_t dim);
    virtual ~WaveletGrid();

    virtual const char *getType();

    virtual GridGenerator *createGridGenerator();

    static Grid *unserialize(std::istream& istr);
};

}
}

#endif /* WAVELETGRID_HPP */
