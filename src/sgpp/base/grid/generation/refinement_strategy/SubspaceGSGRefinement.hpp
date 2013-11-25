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
	typedef std::vector<index_type> ErrorVector;


	SubspaceGSGRefinement(RefinementDecorator* decorator);


	void freeRefineSubspace(GridStorage* storage,RefinementFunctor* functor);

protected:


	SubspaceErrorStorage admissibleSubspaces;
	SubspaceErrorStorage nonAdmissibleSubspaces;
	ErrorVector addedInLastRefinement;

	void filterAdmissibleSubspaces(GridStorage* storage,
								   SubspaceErrorStorage* admissibleSubspaces);

	void updateAdmissibleSubspaces(GridStorage* storage,
								   RefinementFunctor* functor,
								   SubspaceErrorStorage* admissibleSubspaces);

	void hasParentsChecker(GridStorage* storage, index_type& gridPoint, size_t dim, bool result);

	void cleanUpAdmissibleSubspaces(SubspaceErrorStorage* errorStorage,
									ErrorVector* addedInLastStep,
									size_t refinements_num,
									index_type* maxSubspaces
									);

private:


	bool hasParentHelper(GridStorage* storage, index_type& gridPoint, size_t dim);

	void updateAdmissibleSubspacesHelper(GridStorage* storage,
										 RefinementFunctor* functor,
										 SubspaceErrorStorage* errorStorage,
										 ErrorVector* newSubspaces,
										 index_type& storageIndex,
										 size_t dim);
};

} /* namespace base */
} /* namespace sg */
#endif /* SUBSPACEGSGREFINEMENT_HPP_ */
