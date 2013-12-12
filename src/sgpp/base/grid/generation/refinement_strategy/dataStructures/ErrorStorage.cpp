/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "ErrorStorage.hpp"

namespace sg {
namespace base {

ErrorStorage::ErrorStorage(size_t dim):hashErrorStorage(dim) {}

void ErrorStorage::getErrorMap(ErrorMap* errorMap) {
	errorMap = &(this->errorMap);
}

void ErrorStorage::getHashErrorStorage(HashErrorStorage* hashErrorStorage) {
	hashErrorStorage =  &(this->hashErrorStorage);
}

size_t ErrorStorage::push(ErrorType& container) {
	//insert point into hash map
	size_t seq = hashErrorStorage.insert(container);
	//add to multi map sorted by error
	errorMap.insert(std::make_pair(container.getContribPerPoint(),hashErrorStorage[seq]));

	return seq;
}

void ErrorStorage::pop(ErrorType* container) {

	if (errorMap.empty()) {
		//@TODO throw exception
	}
	//get the element with the largest error and copy it to container
	ErrorMap::reverse_iterator mapIter = errorMap.rbegin();
	ErrorType error = *mapIter->second;
	HashErrorStorage::grid_map_iterator hashIter = hashErrorStorage.find(&error);
	container = hashIter->first;

	//delete it from error map and hash storage
	errorMap.erase((++mapIter).base());
	hashErrorStorage.destroy(hashIter->first);

}

void ErrorStorage::peek(ErrorType* container) {

	if (errorMap.empty()) {
		//@TODO throw exception
	}
	//get the element with the largest error and copy it to container
	ErrorMap::reverse_iterator mapIter = errorMap.rbegin();
	HashErrorStorage::grid_map_iterator hashIter = hashErrorStorage.find(mapIter->second);
	container = hashIter->first;
}

void ErrorStorage::popBack(ErrorType* container) {

	if (errorMap.empty()) {
			//@TODO throw exception
		}
		//get the element with the smallest error and copy it to container
		ErrorMap::iterator mapIter = errorMap.begin();
		HashErrorStorage::grid_map_iterator hashIter = hashErrorStorage.find(mapIter->second);
		container = hashIter->first;

		//delete it from error map and hash storage
		errorMap.erase(mapIter);
		hashErrorStorage.destroy(hashIter->first);

}

void ErrorStorage::peekBack(ErrorType* container) {

	if(errorMap.empty()) {
			//@TODO throw exception
	}
		//get the element with the largest error and copy it to container
		ErrorMap::iterator mapIter = errorMap.begin();
		HashErrorStorage::grid_map_iterator hashIter = hashErrorStorage.find(mapIter->second);
		container = hashIter->first;
}

void ErrorStorage::updateErrors(std::vector<size_t>* posVector,
		std::vector<ErrorType> newValVector)
{
	ErrorType tmp;
	for(size_t i = 0; i < posVector->size(); ++i) {

		tmp = *hashErrorStorage[(*posVector)[i]];
		std::pair<ErrorMap::iterator,ErrorMap::iterator> errorObjects = errorMap.equal_range(tmp.getContribPerPoint());

		for(ErrorMap::iterator iter = errorObjects.first; iter != errorObjects.second; ++iter)
		{
			if(iter->second == hashErrorStorage[(*posVector)[i]])
			{
				errorMap.erase(iter);
				hashErrorStorage.destroy(hashErrorStorage[(*posVector)[i]]);
				push(newValVector[i]);
				break;
			}
		}
	}
}


void ErrorStorage::updateErrors(HashErrorStorage* newErrors) {

	for(HashErrorStorage::grid_list_iterator iter = newErrors->begin();
			iter != newErrors->end(); ++iter)
	{
		HashErrorStorage::grid_list_iterator objectToUpdateIter = hashErrorStorage.find(iter->first);

		if(objectToUpdateIter == hashErrorStorage.end())
		{
			//not found insert the error object
			push(*(iter->first));
		}else {
			//found. perform an update
			updateError(*(iter->first));
		}
	}

}

void ErrorStorage::insertErrors(std::vector<ErrorType>* newErrorObjectsVector)
{
	for(std::vector<ErrorType>::iterator iter = newErrorObjectsVector->begin();
			iter != newErrorObjectsVector->end(); ++iter)
	{
		push(*iter);
	}
}

void sg::base::ErrorStorage::insertIntoErrorMap(
	std::vector<size_t>* newErrorObjectsVector) {

	for(std::vector<size_t>::iterator iter = newErrorObjectsVector->begin();
			iter != newErrorObjectsVector->end(); ++iter)
	{
		ErrorType* errorContainer = hashErrorStorage.get(*iter);
		errorMap.insert(std::make_pair(errorContainer->getContribPerPoint(),errorContainer));
	}
}

void sg::base::ErrorStorage::insertAllIntoErrorMap() {

	// go over all indices in the hashErrorStorage
	//and insert all elements into the error map as (contibPerPonint,ErrorContainer) pairs
	for(HashErrorStorage::grid_map_iterator iter= hashErrorStorage.begin();
			iter != hashErrorStorage.end();++iter)
	{
		errorMap.insert(std::make_pair(iter->first->getContribPerPoint(),iter->first));
	}
}


} /* namespace base */
} /* namespace sg */


