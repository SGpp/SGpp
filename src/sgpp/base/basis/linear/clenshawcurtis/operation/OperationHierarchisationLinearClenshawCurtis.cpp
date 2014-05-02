/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/basis/linear/clenshawcurtis/operation/OperationHierarchisationLinearClenshawCurtis.hpp"
#include "base/basis/linear/clenshawcurtis/algorithm_sweep/HierarchisationLinearClenshawCurtis.hpp"
#include "base/basis/linear/clenshawcurtis/algorithm_sweep/DehierarchisationLinearClenshawCurtis.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg
{
namespace base
{

void OperationHierarchisationLinearClenshawCurtis::doHierarchisation(DataVector &node_values)
{
    HierarchisationLinearClenshawCurtis func(storage);
    sweep<HierarchisationLinearClenshawCurtis> s(func, storage);
    
    if (storage->dim() > 1)
    {
        for (size_t i = 0; i < storage->dim(); i++)
        {
            s.sweep1D_Boundary(node_values, node_values, i);
        }
    } else
    {
        s.sweep1D(node_values, node_values, 0);
    }
}

void OperationHierarchisationLinearClenshawCurtis::doDehierarchisation(DataVector &alpha)
{
    DehierarchisationLinearClenshawCurtis func(storage);
    sweep<DehierarchisationLinearClenshawCurtis> s(func, storage);
    
    if (storage->dim() > 1)
    {
        for (size_t i = 0; i < storage->dim(); i++)
        {
            s.sweep1D_Boundary(alpha, alpha, i);
        }
    } else
    {
        s.sweep1D(alpha, alpha, 0);
    }
}

}
}
