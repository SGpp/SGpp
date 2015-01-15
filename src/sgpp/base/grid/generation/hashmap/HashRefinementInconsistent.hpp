/* ****************************************************************************
 * Copyright (C) 2014 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#ifndef HASHREFINEMENTINCONSISTENT_HPP
#define HASHREFINEMENTINCONSISTENT_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/hashmap/HashRefinement.hpp"

namespace sg {
namespace base {

/**
 * Free refinement class for sparse grids
 */
class HashRefinementInconsistent: public HashRefinement {

    using HashRefinement::HashRefinement;

protected:

    /**
     * This method creates a new point on the grid. It checks if some parents or
     * children are needed in other dimensions.
     *
     * @param storage hashmap that stores the gridpoints
     * @param index The point that should be inserted
     */
    void createGridpoint(GridStorage* storage, index_type& index) override;


};
}
}

#endif /* HASHREFINEMENTINCONSISTENT_HPP */
