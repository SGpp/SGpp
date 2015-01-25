/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "SubspaceGSGRefinement.hpp"

namespace sg {
namespace base {

SubspaceGSGRefinement::SubspaceGSGRefinement(AbstractRefinement* refinement, size_t dim):SubspaceRefinement(refinement), availableSubspaces(dim), firstRefinement(true){}

void SubspaceGSGRefinement::freeRefineSubspace(GridStorage* storage,RefinementFunctor* functor)
{

	//nothing there => nothing to refine
	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	//get the HashErrorStorage and a Multimap sorted by error from the ErrorStorage DataStructure.
	HashErrorStorage* subspaceStorage = availableSubspaces.getHashErrorStorage();

	//accumulate error on refinable subspaces
	if (firstRefinement) {

		collectRefinableSubspaces(storage,functor,subspaceStorage);

		//set all refinable subspaces as admissible, where parents in all dimensions exist
		for(HashErrorStorage::grid_map_iterator errorIter = subspaceStorage->begin();
					errorIter != subspaceStorage->end(); ++errorIter)
			{
				ErrorType subspace = *(errorIter->first);
				errorIter->first->setAdmissible(checkAdmissibility(storage,subspace));
			}
		//add all the Subspaces to a map, sorted by error.
		availableSubspaces.insertAllIntoErrorMap();

		firstRefinement = false;

	}else{


		// check if refinement in last step made nonadmissible Subspaces admissible
		for(ErrorVector::iterator lastAddedIter = addedInLastRefinement.begin();
			lastAddedIter !=  addedInLastRefinement.end();
			++lastAddedIter)
		{
			ErrorType neighbour = *lastAddedIter;
			index_t index = 1;
			level_t level = 1;

			for (size_t dim = 0; dim < storage->dim(); ++dim)
			{
				neighbour.get(dim,level,index);
				neighbour.set(dim,level+1,index);

				HashErrorStorage::grid_map_iterator errorIter = subspaceStorage->find(&neighbour);
				if(errorIter!= subspaceStorage->end())
				{
					errorIter->first->setAdmissible(checkAdmissibility(storage,neighbour));
				}
				neighbour = *lastAddedIter;
			}

		}

		updateAdmissibleSubspaces(storage,functor,&addedInLastRefinement,&availableSubspaces);
	}

//	//DEBUG: print all elements
//	std::cout << "\n\navailable subspaces" << "\n";
//	std::cout << "=================================================" << "\n";
//	for(HashErrorStorage::grid_map_iterator errorIter = subspaceStorage->begin(); errorIter != subspaceStorage->end(); ++errorIter)
//	{
//		std::cout << errorIter->first->toString() << std::endl;
//	}


	//-refine all subspaces which satisfy the refinement criteria
	//-empty addedInLastRefinement and insert the subspaces that satisfy the refinement criteria into the storage
	//- remove newly created subspaces from available subspaces.
	refineSubspaceCollection(storage,&availableSubspaces,&addedInLastRefinement,functor);


//	std::cout << "\nremaining subspaces" << "\n";
//		std::cout << "=================================================" << "\n";
//		for(HashErrorStorage::grid_map_iterator errorIter = subspaceStorage->begin(); errorIter != subspaceStorage->end(); ++errorIter)
//		{
//			std::cout << errorIter->first->toString() << std::endl;
//		}
}

void SubspaceGSGRefinement::refineSubspaceCollection(GridStorage* storage,
		ErrorStorage* errorStorage,
		ErrorVector* addedInLastStep,
		RefinementFunctor* functor)
{

		addedInLastStep->clear();

		ErrorMap* errorMap = errorStorage->getErrorMap();

		// now refine all subspaces which satisfy the refinement criteria
		double threshold = functor->getRefinementThreshold();
		size_t refinements_num = functor->getRefinementsNum();
		size_t refined = 0;


		while(refined<refinements_num && !errorMap->empty())
		{
			ErrorType* maxErrorSubspace = errorStorage->peek();
			ErrorType tmp = *maxErrorSubspace;

			if(maxErrorSubspace->isAdmissible() &&
					maxErrorSubspace->getContribPerPoint() > functor->start()
					&& fabs(maxErrorSubspace->getContribPerPoint()) >= threshold)
			{
				createSubspace(storage,tmp,true);
				++refined;
				errorStorage->pop();
				addedInLastStep->push_back(*maxErrorSubspace);
			}

		}

//		//DEBUG: selected
//		std::cout << "\nselected subspaces" << "\n";
//		std::cout << "=================================================" << "\n";
//		for (ErrorVector::iterator iter = addedInLastStep->begin(); iter != addedInLastStep->end(); ++iter)
//		{
//			std::cout << (*iter).toString() << std::endl;
//		}
}

void SubspaceGSGRefinement::updateAdmissibleSubspaces(GridStorage* storage,
												      RefinementFunctor* functor,
												      ErrorVector* addedInLastRefinement,
												      ErrorStorage* availableSubspaces)
{

	HashErrorStorage* errorStorage = availableSubspaces->getHashErrorStorage();
	HashErrorStorage updateSubspaces(storage->dim());

	//go through all subspaces refined in last refinement step.
	for (ErrorVector::iterator addedInLastIter = addedInLastRefinement->begin();
			addedInLastIter != addedInLastRefinement->end(); ++addedInLastIter)
	{
		//create all gridPoint indices in this subspace.
		ErrorType subspace = *addedInLastIter;

		GridPointVector gridPoints;
		ErrorType helper = subspace;
		createAllGridPointsOfSubspace(helper,&gridPoints);
		//here the error contributions will be stored.
		ErrorType errorContribution;


		for(GridPointVector::iterator pointIter = gridPoints.begin(); pointIter != gridPoints.end();++pointIter)
		{
			// functor in each dimension is equal -> error Contribution of each subspace to each neighbouring subspace is equal.
			RefinementFunctor::value_type error =  (*functor)(storage, storage->find(&(*pointIter))->second);

			//add error contrib for left child and for right child
			//careful! += 2*error is wrong, it does not increase the internal contributions counter!
			errorContribution += error;
			errorContribution += error;
		}

		//find all neighbouring subspaces
		index_t index = 1;
		level_t level = 1;

		for(size_t dim = 0; dim < storage->dim(); ++dim)
		{
			ErrorType neighbourSubspace = subspace;
			neighbourSubspace.get(dim,level,index);
			neighbourSubspace.set(dim,level+1,index);

			HashErrorStorage::grid_map_iterator updatedSubspacesIter = updateSubspaces.find(&neighbourSubspace);
			if(updatedSubspacesIter == updateSubspaces.end())
			{
				//not found. insert that subspace.
				neighbourSubspace+=errorContribution;
				updateSubspaces.insert(neighbourSubspace);

			}else{
				//found. update Error Contrib.
				*(updatedSubspacesIter->first)+=errorContribution;
			}
		}
	}

	//update available subspaces
	availableSubspaces->updateErrors(&updateSubspaces);

	//check updated subspaces for admissibility.
	for (HashErrorStorage::grid_map_iterator updatedSubspaceIter = updateSubspaces.begin();
				updatedSubspaceIter != updateSubspaces.end(); ++updatedSubspaceIter)
	{
		HashErrorStorage::grid_map_iterator errorIter = errorStorage->find(updatedSubspaceIter->first);
		errorIter->first->setAdmissible(checkAdmissibility(storage, *(errorIter->first)));
	}
}

} /* namespace base */
} /* namespace sg */


