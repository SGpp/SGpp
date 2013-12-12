/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef ERRORSTORAGE_HPP_
#define ERRORSTORAGE_HPP_

#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/generation/refinement_strategy/dataStructures/ErrorContainer.hpp"
#include "base/grid/storage/hashmap/HashGridStorage.hpp"
#include "base/grid/storage/hashmap/HashGridIterator.hpp"

#include <map>
#include <vector>

namespace sg {
namespace base {

typedef ErrorContainer<AbstractRefinement::level_t,AbstractRefinement::index_t> ErrorType;
typedef HashGridStorage<ErrorType> HashErrorStorage;
typedef std::multimap<RefinementFunctor::value_type,ErrorType*> ErrorMap;

class ErrorStorage {
public:

	ErrorStorage(size_t dim);

	size_t push(ErrorType& container);

	void pop(ErrorType* container);

	void peek(ErrorType* container);

	void popBack(ErrorType* container);

	void peekBack(ErrorType* container);

	void getErrorMap(ErrorMap* errorMap);

	void getHashErrorStorage(HashErrorStorage* hashErrorStorage);

	void updateErrors(std::vector<size_t>* posVector, std::vector<ErrorType> newValVector);

	void updateError(ErrorType newErrorObject);

	void updateErrors(HashErrorStorage* newErrors);

	void insertErrors(std::vector<ErrorType>* newErrorObjectsVector);

	void insertIntoErrorMap(std::vector<size_t>* newErrorObjectsVector);

	void insertAllIntoErrorMap();


protected:

	HashErrorStorage hashErrorStorage;
	ErrorMap errorMap;

};

} /* namespace base */
} /* namespace sg */

#endif /* ERRORSTORAGE_HPP_ */
