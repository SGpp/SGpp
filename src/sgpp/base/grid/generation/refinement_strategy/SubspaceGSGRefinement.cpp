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
	if (admissibleSubspaces.empty()) {
		collectRefinableSubspaces(storage,functor,&admissibleSubspaces);
		filterAdmissibleSubspaces(storage,&admissibleSubspaces);
	}else{
		updateAdmissibleSubspaces(storage,functor,&admissibleSubspaces);
	}

	//select refinements_num highest indicator subspaces
	selectHighestErrorSubspaces(&admissibleSubspaces,refinements_num,maxSubspaces,errorsPerSubspace);
	//delete them from admissible subspaces & add them to subspaces created in the last run.
	cleanUpAdmissibleSubspaces(&admissibleSubspaces,&addedInLastRefinement,refinements_num,maxSubspaces);

	//refine all subspaces which satisfy the refinement criteria
	refineSubspaceCollection(storage,functor,refinements_num,maxSubspaces,errorsPerSubspace);

	delete [] maxSubspaces;
	delete [] errorsPerSubspace;
}

void SubspaceGSGRefinement::filterAdmissibleSubspaces(GridStorage* storage,
		SubspaceErrorStorage* admissibleSubspaces)
{
	std::cout << "filtering admissible subspaces\n";
	for(SubspaceErrorStorage::iterator errorIter = admissibleSubspaces->begin();
			errorIter != admissibleSubspaces->end(); ++errorIter)
	{
		index_type index = errorIter->first;
		std::cout << "checking parents for subspace" << ((index_type) (errorIter->first)).toString() << "\n";
		bool hasParentsInAllDim = true;
		hasParentsChecker(storage,index,0,hasParentsInAllDim);
		errorIter->second.setAdmissible(hasParentsInAllDim);
	}
}

void SubspaceGSGRefinement::updateAdmissibleSubspaces(GridStorage* storage,
		RefinementFunctor* functor, SubspaceErrorStorage* admissibleSubspaces) {

	std::cout << "updating admissible subspaces \n";

	ErrorVector newSubspaces;
	//go through all subspaces refined in last refinement step.
	for (ErrorVector::iterator addedInLastIter = addedInLastRefinement.begin();
			addedInLastIter != addedInLastRefinement.end(); addedInLastIter++) {

		//check children in all dim and either add error to existing, non admissible subspaces or create new ones
		std::cout << "calling update helper on " << ((index_type) (*addedInLastIter)).toString() << "\n";
		updateAdmissibleSubspacesHelper(storage,functor,admissibleSubspaces,&newSubspaces,*addedInLastIter,0);
	}

	//check added subspaces for admissibility.
	std::cout << "checking subspaces for admissibility \n";
	for (ErrorVector::iterator newSubspaceIter = newSubspaces.begin();
				newSubspaceIter != addedInLastRefinement.end(); ++newSubspaceIter) {
		SubspaceErrorStorage::iterator errorIter = admissibleSubspaces->find(*newSubspaceIter);
		index_type subspaceIndex = errorIter->first;
		bool hasParentsInAllDim = true;
		hasParentsChecker(storage,subspaceIndex,0,hasParentsInAllDim);
		std::cout << hasParentsInAllDim;
		std::cout << ((index_type) errorIter->first).toString() << errorIter->second.toString();
		errorIter->second.setAdmissible(hasParentsInAllDim);
		std::cout << "set admissible";
		std::cout << (newSubspaceIter != addedInLastRefinement.end());
		std::cout << addedInLastRefinement.size();
	}
	std::cout << "done";
}

void SubspaceGSGRefinement::hasParentsChecker(GridStorage* storage, index_type& storageIndex, size_t dim,bool result)
{
	//go through every dimension
		cout <<"picking dim " << dim << "<" << storageIndex.dim() << " ?\n";
		if (dim < storageIndex.dim()) {

			//get level of subspace on that dimension
			index_t index;
			level_t level;
			storageIndex.get(dim,level,index);
			index = 1;

			cout << "checking level" << level << " , index " << index << "\n";
			//iterate over all allowed indices on that level, in that dimension
			cout << index << "<" << static_cast <size_t>( 1 << level) << "?\n";
			while(index < static_cast <size_t>( 1 << level)){

				//set gridpoint's index accordingly
				storageIndex.set(dim,level,index);
				std::cout << "currently on  " << dim << ", " << level << ", " << index << "\n";
				//recursive call, so that we iterate over all indices on all levels in all dimensions
				bool test = true;
				hasParentsChecker(storage,storageIndex,dim+1,test);
				result = result && test;
				//move to next admissible
				index=index+2;
				//std::cout << "returning at 1 (" << storageIndex.toString()<< ")\n";
				//return result;
			}

		} else {
			//reached end of recursion scheme. add a new grid point.
			std::cout << "starting has parentHelper for " << storageIndex.toString() << ", on " << dim << "\n";
			//return hasParentHelper(storage,storageIndex,dim);
			result = result && hasParentHelper(storage,storageIndex,dim-1);
		}
}

bool SubspaceGSGRefinement::hasParentHelper(GridStorage* storage, AbstractRefinement::index_type& gridPoint, size_t dim)
{
	AbstractRefinement::index_type parent = gridPoint;
	bool parentOnLevel1 = getParentLevelAndIndex(&parent,dim);
	std::cout << "parent of " << gridPoint.toString() << "in dim " << dim << " is " << parent.toString() << "\n";
	std::cout << "parent on level 1 ? - " << parentOnLevel1 << "\n";
	if(!parentOnLevel1)
	{
		GridStorage::grid_map_iterator parentIter = storage->find(&parent);
		std::cout << gridPoint.toString() << " has parent ?" << (parentIter != storage->end()) << "\n";
		return parentIter != storage->end();

	}else{
		std::cout << "Grid point " << gridPoint.toString() << "has parent on level 1 in dim " << dim << "\n";
		return true;
	}

}

void SubspaceGSGRefinement::cleanUpAdmissibleSubspaces(SubspaceErrorStorage* errorStorage,
													   ErrorVector* addedInLastStep,
													   size_t refinements_num,
													   index_type* maxSubspaces)
{
	std::cout << "celaning lastStep";
	addedInLastStep->clear();
	std::cout << "starting cleanup";
	for (size_t i = 0; i < refinements_num; ++i) {
		errorStorage->erase(maxSubspaces[i]);
		addedInLastStep->push_back(maxSubspaces[i]);
	}

}

void SubspaceGSGRefinement::updateAdmissibleSubspacesHelper(GridStorage* storage,
															RefinementFunctor* functor,
															SubspaceErrorStorage* errorStorage,
															ErrorVector* newSubspaces,
															index_type& storageIndex,
															size_t dim)
{
	//go through every dimension
	cout <<"picking dim " << dim << "<" << storageIndex.dim() << " ?\n";
	if (dim < storageIndex.dim()) {

		//get level of subspace on that dimension
		index_t index;
		level_t level;
		storageIndex.get(dim,level,index);
		index = 1;

		std::cout << "checking level" << level << " , index " << index << "\n";
		//iterate over all allowed indices on that level, in that dimension
		std::cout << index << "<" << static_cast <size_t>( 1 << level) << "?\n";
		while(index < static_cast <size_t>( 1 << level)){

			//set gridpoint's index accordingly
			storageIndex.set(dim,level,index);
			//cout << "currently on  " << dim << ", " << level << ", " << index << "\n";
			//recursive call, so that we iterate over all indices on all levels in all dimensions
			updateAdmissibleSubspacesHelper(storage,functor,errorStorage,newSubspaces,storageIndex,dim+1);
			//move to next admissible
			index=index+2;
		}

	} else {
		//reached end of recursion scheme. add error to the subspace.
		std::cout <<"updating error from point " << storageIndex.toString() << "\n";
		//get the error for gridpoint
		GridStorage::grid_map_iterator gridIter = storage->find(&storageIndex);
		RefinementFunctor::value_type error =  (*functor)(storage, gridIter->second);

		//get its coordinates
		index_t index;
		level_t level;
		storageIndex.get(dim,level,index);

		//for every dimension, get the neighbouring subspace.
		for(size_t d = 0; d<storage->dim();++d)
		{
			//find neighbouring subspace.
			index_type neighbourSubspace = storageIndex;
			storageIndex.get(d,level,index);
			std::cout << "set neighbour subspace to " << neighbourSubspace.toString()<< "\n";
			std::cout << level +1 << std::endl;
			neighbourSubspace.set(d,level+1,1);
			std::cout << "set neighbour subspace to " << neighbourSubspace.toString()<< "\n";
			SubspaceErrorStorage::iterator neighbourIter = errorStorage->find(neighbourSubspace);
			std::cout << "searching for neighbouring subspace " << neighbourSubspace.toString() <<  " in dim " << d <<"\n";
			if(neighbourIter == errorStorage->end())
			{
				std::cout << "not found. creating a new subspace\n";
				//if the neighbouring subspace does not exist, we will create it.
				resetIndexVector(&neighbourSubspace);
				std::pair<SubspaceErrorStorage::iterator,bool> insertionLocation =
						errorStorage->insert(std::make_pair(neighbourSubspace,ErrorContainer()));
				newSubspaces->push_back(neighbourSubspace);
				neighbourIter = insertionLocation.first;
			}

			// add error left child;
			 std::cout << "adding error to " <<  neighbourIter->second.toString() << "\n";
			(neighbourIter->second)+= error;
			// add error right child
			(neighbourIter->second)+= error;
			std::cout << "updated error to " <<  neighbourIter->second.toString() << "\n";
		}
	}
}

} /* namespace base */
} /* namespace sg */
