/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include <algorithm>

#include "opt/basis/linear/clenshawcurtis/operation/OperationMultipleHierarchisationLinearClenshawCurtis.hpp"
#include "opt/basis/linear/clenshawcurtis/operation/OperationEvalLinearClenshawCurtis.hpp"
#include "opt/sle/system/Hierarchisation.hpp"
#include "opt/sle/solver/Auto.hpp"

namespace sg
{
namespace opt
{

void OperationMultipleHierarchisationLinearClenshawCurtis::doHierarchisation(
        std::vector<base::DataVector *> node_values)
{
    sle::system::Hierarchisation system(grid);
    sle::solver::Auto solver;
    std::vector<std::vector<double> > B;
    std::vector<std::vector<double> > X;
    
    for (size_t i = 0; i < node_values.size(); i++)
    {
        B.push_back(std::vector<double>(node_values[i]->getPointer(),
                                        node_values[i]->getPointer() + node_values[i]->getSize()));
    }
    
    if (solver.solve(system, B, X))
    {
        for (size_t i = 0; i < node_values.size(); i++)
        {
            std::copy(X[i].begin(), X[i].begin() + node_values[i]->getSize(),
                      node_values[i]->getPointer());
        }
    }
}

void OperationMultipleHierarchisationLinearClenshawCurtis::doDehierarchisation(
        std::vector<base::DataVector *> alpha)
{
    base::GridStorage *storage = grid.getStorage();
    const size_t d = storage->dim();
    OperationEvalLinearClenshawCurtis op_eval(storage, grid.getCosineTable());
    std::vector<double> node_values_vector(storage->size(), 0.0);
    std::vector<double> x(d, 0.0);
    
    for (size_t i = 0; i < storage->size(); i++)
    {
        for (size_t j = 0; j < storage->size(); j++)
        {
            base::GridIndex *gp = storage->get(j);
            
            for (size_t t = 0; t < d; t++)
            {
                x[t] = gp->abs(t);
            }
            
            node_values_vector[j] = op_eval.eval(*alpha[i], x);
        }
        
        std::copy(node_values_vector.begin(), node_values_vector.begin() + alpha[i]->getSize(),
                  alpha[i]->getPointer());
    }
}

}
}
