// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef HASHREFINEMENTINCONSISTENT_HPP
#define HASHREFINEMENTINCONSISTENT_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Free refinement class for sparse grids
 */
class HashRefinementInconsistent: public HashRefinement {

    //using HashRefinement::HashRefinement;

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
