/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_GRID_WAVELETGRID_HPP
#define SGPP_OPT_GRID_WAVELETGRID_HPP

#include <iostream>

#include "base/grid/Grid.hpp"

namespace sg
{
namespace opt
{

/**
 * Noboundary grid with wavelet basis functions.
 */
class WaveletGrid : public base::Grid
{
public:
    /**
     * Constructor.
     * 
     * @param dim       number of dimensions
     */
    WaveletGrid(size_t dim);
    
    /**
     * Destructor.
     */
    virtual ~WaveletGrid();
    
    /**
     * @return  identifying grid type string
     */
    virtual const char *getType();
    
    /**
     * @return grid generator for this grid type
     */
    virtual base::GridGenerator *createGridGenerator();
    
    /**
     * @param istr  input stream containing the serialization
     * @return      pointer to newly generated deserialized grid
     */
    static base::Grid *unserialize(std::istream &istr);
    
protected:
    /**
     * Deserialization constructor.
     * 
     * @param istr  serialized grid
     */
    WaveletGrid(std::istream &istr);
};

}
}

#endif
