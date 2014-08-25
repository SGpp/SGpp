/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPERATION_OPFACTORY_HPP
#define SGPP_OPT_OPERATION_OPFACTORY_HPP

#include "base/grid/Grid.hpp"

#include "base/operation/OperationEval.hpp"
#include "opt/operation/OperationEvalGradient.hpp"
#include "opt/operation/OperationEvalHessian.hpp"
#include "opt/operation/OperationEvalPartialDerivative.hpp"
#include "opt/operation/OperationMultipleHierarchisation.hpp"

namespace sg
{
namespace opt
{

/**
 * Creates a base::OperationEval for the given sg::opt grid.
 * Don't forget to delete the object after use.
 * 
 * @param grid  sparse grid
 * @return      pointer to a base::OperationEval object for the grid
 */
base::OperationEval *createOperationEval(base::Grid &grid);

/**
 * Creates a OperationEvalGradient for the given sg::opt grid (except for the linear ones).
 * Don't forget to delete the object after use.
 * 
 * @param grid  sparse grid
 * @return      pointer to a OperationEvalGradient object for the grid
 */
OperationEvalGradient *createOperationEvalGradient(base::Grid &grid);

/**
 * Creates a OperationEvalHessian for the given sg::opt grid (except for the linear ones).
 * Don't forget to delete the object after use.
 * 
 * @param grid  sparse grid
 * @return      pointer to a OperationEvalHessian object for the grid
 */
OperationEvalHessian *createOperationEvalHessian(base::Grid &grid);

/**
 * Creates a OperationMultipleHierarchisation for the given sg::opt grid .
 * Don't forget to delete the object after use.
 * 
 * @param grid  sparse grid
 * @return      pointer to a OperationMultipleHierarchisation object for the grid
 */
OperationMultipleHierarchisation *createOperationMultipleHierarchisation(base::Grid &grid);

/**
 * Creates a base::OperationEvalPartialDerivative for the given sg::opt grid
 * (except for the linear ones).
 * Don't forget to delete the object after use.
 * 
 * @param grid  sparse grid
 * @return      pointer to a base::OperationEvalPartialDerivative object for the grid
 */
OperationEvalPartialDerivative *createOperationEvalPartialDerivative(base::Grid &grid);

}
}

#endif
