/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef SUBSPACEGSGREFINEMENT_HPP_
#define SUBSPACEGSGREFINEMENT_HPP_

#include "SubspaceRefinement.hpp"

namespace sg {
namespace base {

class SubspaceGSGRefinement: public SubspaceRefinement {
public:
	typedef std::vector<index_type> SubspaceVector;
	typedef std::vector<index_type> GridPointVector;


	SubspaceGSGRefinement(RefinementDecorator* decorator);


	void freeRefineSubspace(GridStorage* storage,RefinementFunctor* functor);

protected:


	SubspaceErrorStorage availableSubspaces;
	SubspaceVector addedInLastRefinement;

	bool checkAdmissible(GridStorage* storage, index_type& subspace);

	virtual void refineSubspaceCollection(GridStorage* storage,
			SubspaceErrorStorage* errorStorage,
			SubspaceVector* addedInLastStep,
			RefinementFunctor* functor,
			size_t refinements_num,
			index_type* maxErrorSubspaces,
			RefinementFunctor::value_type* maxErrorValues);

	void updateAdmissibleSubspaces(GridStorage* storage,
			RefinementFunctor* functor,
			SubspaceVector* addedInLastRefinement,
			SubspaceErrorStorage* admissibleSubspaces);

	void createAllGridPointsOfSubspace(index_type& subspace, GridPointVector* gridPoints);

private:

	void createAllGridPointsOfSubspaceHelper(GridPointVector* gridPoints,index_type& gridObject,size_t dim);
};

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACEGSGREFINEMENT_HPP_ */
