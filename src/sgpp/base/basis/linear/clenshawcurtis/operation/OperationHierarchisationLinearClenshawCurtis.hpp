/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef OPERATIONHIERARCHISATIONLINEARCLENSHAWCURTIS_HPP
#define OPERATIONHIERARCHISATIONLINEARCLENSHAWCURTIS_HPP

#include "base/operation/OperationHierarchisation.hpp"
#include "base/grid/GridStorage.hpp"

namespace sg
{
namespace base
{

class OperationHierarchisationLinearClenshawCurtis : public OperationHierarchisation
{
public:
    OperationHierarchisationLinearClenshawCurtis(GridStorage *storage) : storage(storage) {}
    virtual ~OperationHierarchisationLinearClenshawCurtis() {}
    
    virtual void doHierarchisation(DataVector &node_values);
    virtual void doDehierarchisation(DataVector &alpha);
    
protected:
    GridStorage* storage;
};

}
}

#endif /* OPERATIONHIERARCHISATIONLINEARCLENSHAWCURTIS_HPP */
