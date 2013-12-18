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

	ErrorStorage();

	size_t push(ErrorType& container);

	ErrorType pop();

	ErrorType popNext();

	ErrorType* peek();

	ErrorType* peekNext();

	ErrorType popBack();

	ErrorType popBackNext();

	ErrorType* peekBack();

	ErrorType* peekBackNext();

	ErrorMap* getErrorMap();

	HashErrorStorage* getHashErrorStorage();

	void addToOldError(ErrorType& newErrorObject);

	void updateErrors(HashErrorStorage* newErrors);

	void insertAllIntoErrorMap();


protected:

	HashErrorStorage hashErrorStorage;
	ErrorMap errorMap;

private:

	ErrorType popSubroutine(ErrorMap::iterator iter);
};

} /* namespace base */
} /* namespace sg */

#endif /* ERRORSTORAGE_HPP_ */
