/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_SLE_SYSTEM_CLONEABLE_HPP
#define SGPP_OPT_SLE_SYSTEM_CLONEABLE_HPP

#include "opt/sle/system/System.hpp"
#include "opt/tools/SmartPointer.hpp"

namespace sg
{
namespace opt
{
namespace sle
{
namespace system
{

/**
 * Abstract class for "cloneable" linear systems.
 * This class is needed in the case that matrix entry lookups are not possible concurrently
 * (e.g. for hierarchisation systems with Clenshaw-Curtis grids).
 */
class Cloneable : public System
{
public:
    /**
     * Constructor.
     */
    Cloneable() : System()
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~Cloneable()
    {
    }
    
    /**
     * Pure virtual method for cloning the linear system.
     * It should return a pointer to the cloned object and it's used for parallel computations
     * (e.g. the getMatrixEntry() method might not be thread-safe).
     * 
     * @return smart pointer to cloned object
     */
    virtual tools::SmartPointer<Cloneable> clone() = 0;
    
    /**
     * @return whether this system derives from Cloneable or not (true)
     */
    bool isCloneable() const
    {
        return true;
    }
};

}
}
}
}

#endif
