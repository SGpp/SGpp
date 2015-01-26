/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef ERRORSTORAGE_HPP_
#define ERRORSTORAGE_HPP_

#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/dataStructures/ErrorContainer.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridIterator.hpp>

#include <map>
#include <vector>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

typedef ErrorContainer<AbstractRefinement::level_t,AbstractRefinement::index_t> ErrorType;
typedef HashGridStorage<ErrorType> HashErrorStorage;
typedef std::multimap<RefinementFunctor::value_type,ErrorType*> ErrorMap;


/**
 * ErrorStorage provides hash map based storage, to store error indicators for grid points or hierarchical subspaces
 * and a tree based structure to sort the error indicators according to the size of their indicators.
 */
class ErrorStorage {
public:

	/**
	 * Constructor
	 *
	 * @param dim dimensionality of the grid used
	 */
	ErrorStorage(size_t dim);

	/**
	 * Default constructor
	 */
	ErrorStorage();

	/**
	 * push
	 * adds an error container to the data structuces in the error storage.
	 *
	 * @param container - the container that shoul be added
	 * @return size_t the index of the container in the hash map.
	 */
	size_t push(ErrorType& container);

	/**
	 * return the error indicator with the highest error and remove it from
	 * the error storage
	 *
	 * @return ErrorType an error container with the highest error indicator.
	 */
	ErrorType pop();

	/**
	 * return the error indicator with the highest error which is admissible and remove it from
	 * the error storage
	 *
	 * @return ErrorType an admissible error container  with the highest error indicator.
	 */
	ErrorType popNext();

	/**
	 * return the error indicator with the highest error without removing it from
	 * the error storage
	 *
	 * @return ErrorType an error container with the highest error indicator.
	 */
	ErrorType* peek();

	/**
	 * return the error indicator with the highest error which is admissible without removing it from
	 * the error storage
	 *
	 * @return ErrorType an admissible error container with the highest error indicator.
	 */
	ErrorType* peekNext();

	/**
	 * return the error indicator with the smallest error and remove it from
	 * the error storage
	 *
	 * @return ErrorType an error container with the smallest error indicator.
	 */
	ErrorType popBack();

	/**
	 * return the error indicator with the smallest error which is admissible and remove it from
	 * the error storage
	 *
	 * @return ErrorType an admissible error container  with the smallest error indicator.
	 */
	ErrorType popBackNext();

	/**
	 * return the error indicator with the smallest error without removing it from
	 * the error storage
	 *
	 * @return ErrorType an error container with the highest error indicator.
	 */
	ErrorType* peekBack();

	/**
	 * return the error indicator with the smallest error which is admissible without removing it from
	 * the error storage
	 *
	 * @return ErrorType an admissible error container with the smallest error indicator.
	 */
	ErrorType* peekBackNext();

	/**
	 * getter for the errorMap
	 *
	 * @returns the underlying multimap, which stores the sorted
	 * error indicators from the the error storage
	 */
	ErrorMap* getErrorMap();

	/**
	 * getter for the hashErrorStorage
	 *
	 * @returns the underlying hash map of the error storage
	 */
	HashErrorStorage* getHashErrorStorage();

	/**
	 * addToOldError searches the ErrorContainer belonging to the same
	 * grid point / subspace and adds its indicator values and contrubutor count
	 * to the already stored indicator for the grid point / subspace
	 *
	 * @param newErrorObject ErrorContainer who's values should be added to an already
	 * existing errorContainer in the error storage.
	 */
	void addToOldError(ErrorType& newErrorObject);

	/**
	 * updateErrors inserts all ErrorContainers from the argument into the existing error storage.
	 * If ErrorContainers in from the error storage already exist, the values from the argument are added to the existing values.
	 *
	 * @param newErrors HashErrorStorage containing ErrorContainers to be added to the errorStorage.
	 */
	void updateErrors(HashErrorStorage* newErrors);

	/**
	 * Inserts all ErrorContainers from the error storage's internal
	 * hash map into the error storage's internal map that is sorted by error indicators
	 */
	void insertAllIntoErrorMap();

	/**
	 * rebuilds the error storage's internal map that is sorted by error indicators
	 */
	void rebuildErrorMap();


protected:

	//hash map based storage for error indicators
	HashErrorStorage hashErrorStorage;
	//tree based storage for error indicators sorted by error indicators, highest error indicators are first.
	ErrorMap errorMap;

private:

	/**
	 * Subroutine shared by pop methods.
	 *
	 * @param iter ErrorMap::iterator to the element from the error map, that should be returned and removed
	 * @return the ErrorContainer at the iter position.
	 */
	ErrorType popSubroutine(ErrorMap::iterator iter);
};

} /* namespace base */
} /* namespace SGPP */

#endif /* ERRORSTORAGE_HPP_ */
