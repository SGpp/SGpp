/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef SGPP_OPT_OPERATION_OPFACTORY_HPP
#define SGPP_OPT_OPERATION_OPFACTORY_HPP

#include "base/grid/Grid.hpp"
#include "opt/operation/OperationMultipleHierarchisation.hpp"

namespace sg {
  namespace op_factory {

    /**
     * Creates a OperationMultipleHierarchisation for the given sg::opt grid .
     * Don't forget to delete the object after use.
     *
     * @param grid  sparse grid
     * @return      pointer to a OperationMultipleHierarchisation object for the grid
     */
    opt::OperationMultipleHierarchisation* createOperationMultipleHierarchisation(base::Grid& grid);

  }
}

#endif
