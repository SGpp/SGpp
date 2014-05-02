/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef WAVELETBOUNDARYGRID_HPP
#define WAVELETBOUNDARYGRID_HPP

#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg
{
namespace base
{

class WaveletBoundaryGrid : public Grid
{
protected:
    WaveletBoundaryGrid(std::istream& istr);
    
public:
    WaveletBoundaryGrid(size_t dim);
    virtual ~WaveletBoundaryGrid();

    virtual const char *getType();

    virtual GridGenerator *createGridGenerator();

    static Grid *unserialize(std::istream& istr);
};

}
}

#endif /* WAVELETBOUNDARYGRID_HPP */
