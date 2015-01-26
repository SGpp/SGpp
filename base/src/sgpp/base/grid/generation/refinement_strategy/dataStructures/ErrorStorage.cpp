// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include "ErrorStorage.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

ErrorStorage::ErrorStorage(size_t dim):hashErrorStorage(dim) {}

ErrorStorage::ErrorStorage():hashErrorStorage(0){}

ErrorMap* ErrorStorage::getErrorMap() {
	return &(this->errorMap);
}

HashErrorStorage* ErrorStorage::getHashErrorStorage() {
	return &(this->hashErrorStorage);
}

size_t ErrorStorage::push(ErrorType& container) {
	//insert point into hash map
	size_t seq = hashErrorStorage.insert(container);
	//add to multi map sorted by error
	errorMap.insert(std::make_pair(container.getContribPerPoint(),hashErrorStorage[seq]));

	return seq;
}

ErrorType ErrorStorage::pop() {

	//get the element with the largest error and copy it to container
	for(ErrorMap::reverse_iterator mapIter = errorMap.rbegin();
			mapIter!=errorMap.rend();++mapIter)
	{
		if(mapIter->second->isAdmissible()){
			return popSubroutine(((++mapIter).base()));
		}
	}

	//no admissible available - return just the next subspace.
	return popBackNext();
}

ErrorType ErrorStorage::popNext() {
	return popSubroutine(((++errorMap.rbegin()).base()));
}

ErrorType* ErrorStorage::peek() {

	//get the element with the largest error and copy it to container
	for(ErrorMap::reverse_iterator mapIter = errorMap.rbegin();
			mapIter!=errorMap.rend();++mapIter)
	{
		if(mapIter->second->isAdmissible()){
			return ((++mapIter).base())->second;
		}
	}

	//no admissible available - return just the next subspace.
	return peekNext();
}

ErrorType* ErrorStorage::peekNext() {
	return (++errorMap.rbegin()).base()->second;
}

ErrorType ErrorStorage::popBack()
{
	//get the element with the largest error and copy it to container
	for(ErrorMap::iterator mapIter = errorMap.begin();
			mapIter!=errorMap.end();++mapIter)
	{
		if(mapIter->second->isAdmissible()){
			return popSubroutine(mapIter);
		}
	}

	//no admissible available - return just the next subspace.
	return popBackNext();
}

ErrorType ErrorStorage::popBackNext() {
	ErrorMap::iterator iter = errorMap.begin();
	return popSubroutine(iter);
}

ErrorType* ErrorStorage::peekBack()
{
	//get the element with the largest error and copy it to container
	for(ErrorMap::iterator mapIter = errorMap.begin();
			mapIter!=errorMap.end();++mapIter)
	{
		if(mapIter->second->isAdmissible()){
			return mapIter->second;
		}
	}

	//no admissible available - return just the next subspace.
	return peekBackNext();
}

ErrorType* ErrorStorage::peekBackNext() {
	return errorMap.begin()->second;
}


void ErrorStorage::addToOldError(ErrorType& newErrorObject)
{

	std::cout << "adding to old error " << newErrorObject.toString() << std::endl;

	//find the old errorObject in the storage to get its current error to
	//find the current error in the error map
	//as this is a multimap - multiple error objects could have the same error -> find the correct one.
	HashErrorStorage::grid_map_iterator oldErrorObjIter = hashErrorStorage.find(&newErrorObject);
	ErrorType oldError = *(oldErrorObjIter->first);

	std::pair<ErrorMap::iterator,ErrorMap::iterator> errorRange =
			errorMap.equal_range(oldErrorObjIter->first->getContribPerPoint());

	for (ErrorMap::iterator iter = errorRange.first; iter!=errorRange.second; ++iter)
	{
		if(newErrorObject.equals(*(iter->second)))
		{
			std::cout << "erasing from map: " << (iter->second)->toString();
			errorMap.erase(iter);
			std::cout << " - done" << std::endl;
			break;
		}
	}

	//now update the value in the has map and reinsert into the errorMap.
	size_t seq = oldErrorObjIter->second;
	std::cout << oldErrorObjIter->first->toString() << std::endl;
	newErrorObject+=oldError;
	std::cout << "updating tp" << newErrorObject.toString()<< std::endl;
	hashErrorStorage.update(newErrorObject,seq);

	errorMap.insert(std::make_pair(hashErrorStorage[seq]->getContribPerPoint(),hashErrorStorage[seq]));

}

void ErrorStorage::updateErrors(HashErrorStorage* newErrors) {

	std::cout << "updating errors\n";

	for(HashErrorStorage::grid_map_iterator iter = newErrors->begin();iter != newErrors->end(); ++iter)
	{
		std::cout << "searching for " << iter->first->toString();
		HashErrorStorage::grid_map_iterator objectToUpdateIter = hashErrorStorage.find(iter->first);

		if(objectToUpdateIter == hashErrorStorage.end())
		{
			std::cout << "- not found. adding new.\n";
			//not found insert the error object
			push(*(iter->first));
		}else {
			//found. perform an update
			std::cout << "- found. updating. \n";
			addToOldError(*(iter->first));
		}
	}
}


void ErrorStorage::insertAllIntoErrorMap() {

	// go over all indices in the hashErrorStorage
	//and insert all elements into the error map as (contibPerPonint,ErrorContainer) pairs
	for(HashErrorStorage::grid_map_iterator iter= hashErrorStorage.begin();
			iter != hashErrorStorage.end();++iter)
	{
		errorMap.insert(std::make_pair(iter->first->getContribPerPoint(),iter->first));
	}
}

void ErrorStorage::rebuildErrorMap() {

	//empty the error map and rebuild it from scratch
	errorMap.clear();
	insertAllIntoErrorMap();
}

ErrorType ErrorStorage::popSubroutine(ErrorMap::iterator iter) {
	if (errorMap.empty()) {
			//@TODO throw exception
		}
		//get the element with the smallest error and copy it to container
		HashErrorStorage::grid_map_iterator hashIter = hashErrorStorage.find(iter->second);
		ErrorType container = *(hashIter->first);

		//delete it from error map and hash storage
		errorMap.erase(iter);
		std::list<size_t> itemsToRemove;
		itemsToRemove.push_back(hashIter->second);
		hashErrorStorage.deletePoints(itemsToRemove);

		return container;
}



} /* namespace base */
} /* namespace SGPP */

