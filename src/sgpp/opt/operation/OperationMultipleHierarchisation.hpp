/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPERATION_OPERATIONMULTIPLEHIERARCHISATION_HPP
#define SGPP_OPT_OPERATION_OPERATIONMULTIPLEHIERARCHISATION_HPP

#include <vector>

#include "base/datatypes/DataVector.hpp"
#include "base/operation/OperationHierarchisation.hpp"

namespace sg
{
namespace opt
{

/**
 * Abstract operation for hierarchisation and dehierarchisation for multiple sets
 * of function values at the grid nodes.
 */
class OperationMultipleHierarchisation : public base::OperationHierarchisation
{
public:
    /**
     * Constructor.
     */
    OperationMultipleHierarchisation()
    {
    }
    
    /**
     * Virtual destructor.
     */
    virtual ~OperationMultipleHierarchisation()
    {
    }
    
    /**
     * Virtual method for hierarchising for one set of function values.
     * Defaults to calling the other doHierarchisation() by doing two copy operations.
     * 
     * @param[in,out] node_values   before: vector of function values at the grid points,
     *                              after: vector of hierarchical coefficients
     */
    virtual void doHierarchisation(base::DataVector &node_values)
    {
        std::vector<base::DataVector *> node_values_vec;
        node_values_vec.push_back(&node_values);
        doHierarchisation(node_values_vec);
    }
    
    /**
     * Virtual method for dehierarchising for one set of function values.
     * Defaults to calling the other doDehierarchisation() by doing two copy operations.
     * 
     * @param[in,out] node_values   before: vector of hierarchical coefficients,
     *                              after: vector of function values at the grid points
     */
    virtual void doDehierarchisation(base::DataVector &alpha)
    {
        std::vector<base::DataVector *> alpha_vec;
        alpha_vec.push_back(&alpha);
        doDehierarchisation(alpha_vec);
    }
    
    /**
     * Pure virtual method for hierarchising for multiple sets of function values.
     * 
     * @param[in,out] node_values   before: vector of function values at the grid points,
     *                              after: vector of hierarchical coefficients
     */
    virtual void doHierarchisation(std::vector<base::DataVector *> node_values) = 0;
    
    /**
     * Pure virtual method for dehierarchising for multiple sets of coefficients.
     * 
     * @param[in,out] alpha         before: vector of hierarchical coefficients,
     *                              after: vector of function values at the grid points
     */
    virtual void doDehierarchisation(std::vector<base::DataVector *> alpha) = 0;
};

}
}

#endif
