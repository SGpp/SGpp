/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_GRID_BSPLINECLENSHAWCURTISGRID_HPP
#define SGPP_OPT_GRID_BSPLINECLENSHAWCURTISGRID_HPP

#include <iostream>

#include "base/grid/Grid.hpp"
#include "opt/tools/CosineTable.hpp"

namespace sg
{
namespace opt
{

/**
 * Clenshaw-Curtis grid with B-spline basis functions.
 */
class BsplineClenshawCurtisGrid : public base::Grid
{
public:
    /**
     * Constructor.
     * 
     * @param dim           number of dimensions
     * @param degree        B-spline degree
     * @param cosine_table  lookup cosine table (optional)
     */
    BsplineClenshawCurtisGrid(size_t dim, size_t degree,
        const tools::CosineTable *cosine_table = NULL);
    
    /**
     * Destructor.
     */
    virtual ~BsplineClenshawCurtisGrid();
    
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
    
    /**
     * @return  cosine table
     */
    virtual const tools::CosineTable *getCosineTable() const;
    
    /**
     * @param cosine_table  new cosine table
     */
    virtual void setCosineTable(const tools::CosineTable *cosine_table);
    
protected:
    /// B-spline degree
    size_t degree;
    /// cosine table
    const tools::CosineTable *cosine_table;
    
    /**
     * Deserialization constructor.
     * 
     * @param istr  serialized grid
     */
    BsplineClenshawCurtisGrid(std::istream &istr);
};

}
}

#endif
