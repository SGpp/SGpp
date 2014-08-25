/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_GRID_MODIFIEDBSPLINEGRID_HPP
#define SGPP_OPT_GRID_MODIFIEDBSPLINEGRID_HPP

#include <iostream>

#include "base/grid/Grid.hpp"

namespace sg
{
namespace opt
{

/**
 * Noboundary grid with modified B-spline basis functions.
 */
class ModBsplineGrid : public base::Grid
{
public:
    /**
     * Constructor.
     * 
     * @param dim       number of dimensions
     * @param degree    B-spline degree
     */
    ModBsplineGrid(size_t dim, size_t degree);
    
    /**
     * Destructor.
     */
    virtual ~ModBsplineGrid();
    
    /**
     * @return  identifying grid type string
     */
    virtual const char *getType();
    
    /**
     * @return grid generator for this grid type
     */
    virtual base::GridGenerator *createGridGenerator();
    
    /**
     * @param[out] ostr     output stream as target of serialization
     */
    virtual void serialize(std::ostream &ostr);
    
    /**
     * @param istr  input stream containing the serialization
     * @return      pointer to newly generated deserialized grid
     */
    static base::Grid *unserialize(std::istream &istr);
    
    /**
     * @return  B-spline degree
     */
    virtual size_t getDegree();
    
protected:
    /// B-spline degree
    size_t degree;
    
    /**
     * Deserialization constructor.
     * 
     * @param istr  serialized grid
     */
    ModBsplineGrid(std::istream &istr);
};

}
}

#endif
