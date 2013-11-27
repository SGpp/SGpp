/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "SubspaceGSGRefinement.hpp"


using namespace std;

namespace sg {
namespace base {

SubspaceGSGRefinement::SubspaceGSGRefinement(RefinementDecorator* decorator):SubspaceRefinement(decorator){}

//@TODO fix error -> subspaces that are below threshold are removed but not refined.
void SubspaceGSGRefinement::freeRefineSubspace(GridStorage* storage,RefinementFunctor* functor)
{

	//nothing there => nothing to refine
	if (storage->size() == 0) {
		throw generation_exception("storage empty");
	}

	// the functor->getRefinementsNum() largest subspaces should be refined.
	size_t refinements_num = functor->getRefinementsNum();
	//subspaces
	index_type* maxSubspaces = new index_type[refinements_num];
	//subspace errors
	RefinementFunctor::value_type* errorsPerSubspace = new RefinementFunctor::value_type[refinements_num];

	//accumulate error on refinable subspaces
	if (availableSubspaces.empty()) {
		collectRefinableSubspaces(storage,functor,&availableSubspaces);

		//set all refinable subspaces as admissible, where parents in all dimensions exist
		for(SubspaceErrorStorage::iterator errorIter = availableSubspaces.begin();
					errorIter != availableSubspaces.end(); ++errorIter)
			{
				index_type subspace = errorIter->first;
				errorIter->second.setAdmissible(checkAdmissible(storage,subspace));
			}

	}else{

		// check if refinement in last step made nonadmissible Subspaces admissible
		for(SubspaceVector::iterator lastAddedIter = addedInLastRefinement.begin();
			lastAddedIter !=  addedInLastRefinement.end();
			++lastAddedIter)
		{
			index_type neighbour = *lastAddedIter;
			index_t index = 1;
			level_t level = 1;

			for (size_t dim = 0; dim < storage->dim(); ++dim)
			{
				neighbour.get(dim,level,index);
				neighbour.set(dim,level+1,index);
				SubspaceErrorStorage::iterator errorIter = availableSubspaces.find(neighbour);
				errorIter->second.setAdmissible(checkAdmissible(storage,neighbour));
				neighbour = *lastAddedIter;
			}
		}

		updateAdmissibleSubspaces(storage,functor,&addedInLastRefinement,&availableSubspaces);
	}

	//DEBUG: print all elements

	std::cout << "\n\navailable subspaces" << "\n";
	std::cout << "=================================================" << "\n";
	for(SubspaceErrorStorage::iterator errorIter = availableSubspaces.begin(); errorIter != availableSubspaces.end(); ++errorIter)
	{
		std::cout << ((index_type) errorIter->first).toString() << " , " << errorIter->second.toString() << "\n";
	}

	//select refinements_num highest indicator subspaces
	selectHighestErrorSubspaces(&availableSubspaces,refinements_num,maxSubspaces,errorsPerSubspace);

	//DEBUG: selected
	std::cout << "\nselected subspaces" << "\n";
	std::cout << "=================================================" << "\n";
	for (size_t i = 0; i < refinements_num; ++i)
	{
		std::cout << ((index_type) maxSubspaces[i]).toString() << " , " << errorsPerSubspace[i] << "\n";
	}

	//-refine all subspaces which satisfy the refinement criteria
	//-empty addedInLastRefinement and insert the subspaces that satisfy the refinement criteria into the storage
	//- remove newly created subspaces from available subspaces.
	refineSubspaceCollection(storage,&availableSubspaces,&addedInLastRefinement,functor,refinements_num,maxSubspaces,errorsPerSubspace);

	delete [] maxSubspaces;
	delete [] errorsPerSubspace;
}

void SubspaceGSGRefinement::refineSubspaceCollection(GridStorage* storage,
		SubspaceErrorStorage* errorStorage,
		SubspaceVector* addedInLastStep,
		RefinementFunctor* functor,
		size_t refinements_num,
		index_type* maxErrorSubspaces,
		RefinementFunctor::value_type* maxErrorValues)
{
	//cout << "refining subspace collection\n";
		RefinementFunctor::value_type maxErrorValue;
		index_type maxErrorSubspaceIndex;
		addedInLastStep->clear();
		// now refine all grid points which satisfy the refinement criteria
		double threshold = functor->getRefinementThreshold();

		for (size_t i = 0; i < refinements_num; i++) {
			maxErrorValue = maxErrorValues[i];
			maxErrorSubspaceIndex = maxErrorSubspaces[i];

			if (maxErrorValue > functor->start() && fabs(maxErrorValue) >= threshold) {
				createSubspace(storage,maxErrorSubspaceIndex);

				//insert the subspaces that satisfy the refinement criteria into the storage
				// and remove newly created subspaces from available subspaces.
				errorStorage->erase(maxErrorSubspaceIndex);
				addedInLastStep->push_back(maxErrorSubspaceIndex);
			}
		}
}

void SubspaceGSGRefinement::updateAdmissibleSubspaces(GridStorage* storage,
												      RefinementFunctor* functor,
												      SubspaceVector* addedInLastRefinement,
												      SubspaceErrorStorage* availableSubspaces)
{

//	std::cout << "updating admissible subspaces \n";

	SubspaceVector newSubspaces;
	//go through all subspaces refined in last refinement step.
	for (SubspaceVector::iterator addedInLastIter = addedInLastRefinement->begin();
			addedInLastIter != addedInLastRefinement->end(); ++addedInLastIter)
	{

		//create all gridPoints in this subspace.
		index_type subspace = *addedInLastIter;
		GridPointVector gridPoints;
		createAllGridPointsOfSubspace(subspace,&gridPoints);

		//for each grid point, go through all dimensions
		for(GridPointVector::iterator pointIter = gridPoints.begin(); pointIter != gridPoints.end();++pointIter)
		{
			// functor for each child is equal. get it.
			RefinementFunctor::value_type error =  (*functor)(storage, storage->find(&(*pointIter))->second);

			//find out if children exist in
			// 1) a refinable, non admissible subspace,
			// 2) a non existing subspace

			index_t index = 1;
			level_t level = 1;

			for(size_t dim = 0; dim < storage->dim(); ++dim)
			{
				//find neighbouring subspace.
				index_type neighbourSubspace = *pointIter;
				(*pointIter).get(dim,level,index);

				neighbourSubspace.set(dim,level+1,1);

				SubspaceErrorStorage::iterator neighbourIter = availableSubspaces->find(neighbourSubspace);

				if(neighbourIter == availableSubspaces->end())
				{

//					std::cout << "did not find " << neighbourSubspace.toString() << "\n";
					//if the neighbouring subspace does not exist, we will create it.
					resetIndexVector(&neighbourSubspace);
					std::pair<SubspaceErrorStorage::iterator,bool> insertionLocation =
							availableSubspaces->insert(std::make_pair(neighbourSubspace,ErrorContainer()));
					newSubspaces.push_back(neighbourSubspace);
					neighbourIter = insertionLocation.first;
				}
//				else
//				{
//					std::cout << "foud " << neighbourSubspace.toString() << "\n";
//				}

				// add error left child;
				(neighbourIter->second)+= error;
				// add error right child
				(neighbourIter->second)+= error;
				//std::cout << "updated error to " <<  neighbourIter->second.toString() << "\n";
			}
		}
	}

	//check added subspaces for admissibility.
	std::cout << "checking subspaces for admissibility \n";
	for (SubspaceVector::iterator newSubspaceIter = newSubspaces.begin();
				newSubspaceIter != newSubspaces.end(); ++newSubspaceIter)
	{
		SubspaceErrorStorage::iterator errorIter = availableSubspaces->find(*newSubspaceIter);
		index_type subspace = errorIter->first;
		errorIter->second.setAdmissible(checkAdmissible(storage, subspace));
		std::cout << ((index_type) errorIter->first).toString() << " ; " <<  errorIter->second.toString() << "\n";
	}


	// DEBUG:: print all new subspaces
}

void SubspaceGSGRefinement::createAllGridPointsOfSubspace(index_type& subspace,
		GridPointVector* gridPoints) {

	createAllGridPointsOfSubspaceHelper(gridPoints,subspace,0);
}

bool SubspaceGSGRefinement::checkAdmissible(GridStorage* storage,index_type& subspace)
{

	index_type gridPoint = subspace;
	bool isAdmissible = true;
	size_t dim = 0;
	//go through all dimensions and check if all parents are availabe
	while(dim < gridPoint.dim() && isAdmissible)
	{
		//get the parent index
		bool childOnLevelOne = getParentLevelAndIndex(&gridPoint,dim);
		//if the child can have a parent index (= the index in dim is not on lvl 1)
		if(!childOnLevelOne)
		{
			//if we can not find the parent in the grid, the subspace is not admissible;
			isAdmissible = (storage->find(&gridPoint) != storage->end());
		}

		gridPoint = subspace;
		++dim;
	}

	return isAdmissible;
}

void SubspaceGSGRefinement::createAllGridPointsOfSubspaceHelper(
		GridPointVector* gridPoints, index_type& storageIndex, size_t dim)
{

	//go through every dimension
	//cout << dim << "<" << storageIndex.dim() << " ?\n";
	if (dim < storageIndex.dim()) {

		//get level of subspace on that dimension
		index_t index;
		level_t level;
		storageIndex.get(dim,level,index);
		index = 1;

		//cout << level << " , " << index << "\n";
		//iterate over all allowed indices on that level, in that dimension
		//cout << index << "<" << static_cast <size_t>( 1 << level) << "?\n";
		while(index < static_cast <size_t>( 1 << level)){

			//set gridpoint's index accordingly
			storageIndex.set(dim,level,index);
			//cout << "currently on  " << dim << ", " << level << ", " << index << "\n";
			//recursive call, so that we iterate over all indices on all levels in all dimensions
			createAllGridPointsOfSubspaceHelper(gridPoints,storageIndex,dim+1);
			//move to next admissible
			index=index+2;
		}

	} else {
		//reached end of recursion scheme. add a new grid point.
//		cout << "vector: adding gridpoint " << storageIndex.toString() << "\n";
		gridPoints->push_back(storageIndex);
	}
}



} /* namespace base */
} /* namespace sg */


