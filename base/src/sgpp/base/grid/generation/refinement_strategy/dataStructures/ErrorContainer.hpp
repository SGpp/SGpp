// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef ErrorContainer_HPP_
#define ErrorContainer_HPP_

#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridIndex.hpp>
#include <sgpp/base/grid/storage/hashmap/HashGridStorage.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {


/**
 * ErrorContainer.
 *
 * Class derived from HashGridIndex, with the additional ability to store information from
 * a refinement indicator.
 * For each grid object we can store
 * -error
 * -amount of contributors that make up the error
 * -admissibility of the gridpoint /subspace
 */
template<class LT,class IT>
class ErrorContainer: public HashGridIndex<LT,IT> {
public:

	/**
	 * Constructor
	 *
	 * initialize error container, and pass all attributes manually.
	 *
	 * @param error the error indicator value for the grid point / subspace
	 * @param contributionCounter tells how many values make up the indicator value.
	 * Utilized for the caluclation of a mean error for subspaces.
	 * @param index the grid point/subspace for which to store the error indicator.
	 *
	 */
	ErrorContainer(RefinementFunctor::value_type error, size_t contributionCounter, HashGridIndex<LT,IT>& index):HashGridIndex<LT,IT>(index)
	{
		this->error = error;
		this->contributionCounter = contributionCounter;
		admissible = true;
	}

	/**
	 * Constructor
	 *
	 * initialize error container, specifying the error for a grid point / subspace.
	 *
	 * @param error the error indicator value for the grid point / subspace
	 * @param index the grid point/subspace for which to store the error indicator
	 */
	ErrorContainer(RefinementFunctor::value_type error, HashGridIndex<LT,IT>& index):HashGridIndex<LT,IT>(index)
			{
		this->error = error;
		this->contributionCounter = 1;
		admissible = true;
			}


	/**
	 * Constructor
	 *
	 * initialize error container for a grid point / subspace
	 *
	 * @param error the error indicator value for the grid point / subspace
	 * @param index the grid point/subspace for which to store the error indicator
	 */
	ErrorContainer(HashGridIndex<LT,IT>& index):HashGridIndex<LT,IT>(index)
			{
		error = 0;
		contributionCounter = 0;
		admissible = true;
			}


	/**
	 * CopyConstructor using reference
	 *
	 * @oaram container who's values to copy to this object.
	 */
	ErrorContainer(const ErrorContainer<LT,IT>& container):HashGridIndex<LT,IT>(container)
			{
		error = container.error;
		contributionCounter = container.contributionCounter;
		admissible = container.admissible;
			}

	/**
	 * True CopyConstructor using pointer
	 *
	 * @oaram container who's values to copy to this object.
	 */
	ErrorContainer(ErrorContainer<LT,IT>* container):HashGridIndex<LT,IT>(container)
			{
		error = container->error;
		contributionCounter = container->contributionCounter;
		admissible = container->admissible;
			}

	/**
	 * Default constructor for empty ErrorContainer.
	 *
	 * Calls default constructor for HashGridIndex and sets all
	 * values to 0 and marks the HashGridIndex as admissible.
	 */
	ErrorContainer():HashGridIndex<LT,IT>()
			{
		error = 0;
		contributionCounter = 0;
		admissible = true;
			}

	/**
	 * Getter method for the error indicator value stored.
	 *
	 * @return error the error indicator value
	 */
	RefinementFunctor::value_type getError()
	{
		return error;
	}


	/**
	 * Getter method for the amount of contributors that make up the indicator value.
	 *
	 * @return contributionCounter the amount of contributors that make up the indicator value.
	 */
	size_t getContributionCounter()
	{
		return contributionCounter;
	}

	/**
	 * Overloaded = operator. allows to set the stored error indicator value to the value
	 * provided. The contributor count is reset to 1.
	 *
	 * @param newError new value to set the error to.
	 */
	void operator= (RefinementFunctor::value_type newError)
	{
		error = newError;
		contributionCounter = 1;
	}


	/**
	 * Overloaded = operator, to support correct deep copys of the parameters
	 * in the error container.
	 *
	 * @param container ErrorContainer containing the values
	 * to set the member variables of the current object to
	 * @return object - a deep copy of the argument object.
	 */
	ErrorContainer<LT,IT>& operator= (ErrorContainer<LT,IT>& container)
	{
		HashGridIndex<LT,IT>::operator =(container);
		error = container.error;
		contributionCounter = container.contributionCounter;
		admissible = container.admissible;

		return *this;
	}

	/**
	 * Overloaded + operator, allows addition of new contributor to error indicator.
	 * The contributor count is automatically increased.
	 *
	 * @param newError
	 * @return ErrorContainer instance with summed error indicators
	 */
	ErrorContainer<LT,IT> operator+ (RefinementFunctor::value_type newError)
	{
		ErrorContainer result = *this;
		result.error += newError;
		result.contributionCounter++;

		return result;
	}

	/**
	 * Overloaded + operator, allows addition of two error Indicators
	 * The error indicator and the contribution counter of the container from the argument
	 * are added to the corresponding variables in the current indicator.
	 *
	 *@param container the container, who's error indicator and
	 *counter to add to the current ErrorContainer instance
	 *@return ErrorContainer instance with summed error indicators
	 */
	ErrorContainer<LT,IT> operator+ (ErrorContainer<LT,IT>& container)
	{
		ErrorContainer result = *this;
		result.error += container.error;
		result.contributionCounter += container.contributionCounter;

		return result;
	}


	/**
	 * Overloaded =+ operator, allows addition of new contributor to error indicator.
	 * The contributor count is automatically increased.
	 *
	 * @param newError indicator value to add
	 * @return an ErrorContainer instance with summed error indicators
	 */
	ErrorContainer<LT,IT> operator+= (RefinementFunctor::value_type newError)
			{
		this->error += newError;
		this->contributionCounter ++;
		return *this;
			}

	/**
	 * Overloaded += operator, allows addition of two error Indicators
	 * The error indicator and the contribution counter of the container from the argument
	 * are added to the corresponding variables in the current indicator.
	 *
	 *@param container the container, who's error indicator and
	 *counter to add to the current ErrorContainer instance
	 *@return ErrorContainer instance with summed error indicators
	 */
	ErrorContainer<LT,IT>& operator+= (ErrorContainer<LT,IT>& container)
			{
		this->error += container.error;
		this->contributionCounter += container.contributionCounter;
		return *this;
			}

	/**
	 * Overloaded - operator, allows subtraction of new contributor to error indicator.
	 * The contributor count is automatically increased.
	 *
	 * @param newError indicator value to subtract
	 * @return an ErrorContainer instance with the input parameter subtracted
	 * from the stored error indicator
	 */
	ErrorContainer<LT,IT> operator- (RefinementFunctor::value_type newError)
	{
		ErrorContainer result = *this;
		result.error -= newError;
		result.contributionCounter++;

		return result;
	}

	/**
	 * Overloaded - operator, allows subtraction of two error Indicators
	 * The error indicator and the contribution counter of the container from the argument
	 * are subtracted from to the corresponding variables in the current indicator.
	 *
	 *
	 *
	 */
	ErrorContainer<LT,IT> operator- (ErrorContainer<LT,IT>& container)
	{
		ErrorContainer result = *this;
		result.error -= container.error;
		result.contributionCounter += container.contributionCounter;

		return result;
	}

	/**
	 * Overloaded == operator, compares error/numContributors for equality with the specified value.
	 *
	 * @param otherValue value to compare to
	 * @return true if error/numContributors is equal to the input parameter
	 */
	bool operator== (RefinementFunctor::value_type otherValue)
			{
		return getContribPerPoint() == otherValue;
			}

	/**
	 * Overloaded == operator, compares error/numContributors for equality with the argument's error/numContributors
	 *
	 * @param container container to compare to
	 * @return true if both error/numContributors are equal
	 */
	bool operator== (ErrorContainer<LT,IT>& container)
			{
		return getContribPerPoint() == container.getContribPerPoint();
			}


	/**
	 * Overloaded > operator, compares if error/numContributors > argument
	 *
	 * @param otherValue value to compare to
	 * @return true if error/numContributors > argument
	 */
	bool operator> (RefinementFunctor::value_type otherValue)
	{
		return getContribPerPoint() > otherValue;
	}


	/**
	 * Overloaded > operator, compares if this container's error/numContributors > argument's error/numContributors
	 *
	 * @param container ErrorContainer to compare to
	 * @return true if container's error/numContributors > argument's error/numContributors
	 */
	bool operator> (const ErrorContainer<LT,IT>& container) const
	{
		return getContribPerPoint() > container.getContribPerPoint();
	}


	/**
	 * Overloaded < operator, compares if error/numContributors < argument
	 *
	 * @param otherValue value to compare to
	 * @return true if container's error/numContributors < argument
	 */
	bool operator< (RefinementFunctor::value_type otherValue)
	{
		return getContribPerPoint() < otherValue;
	}


	/**
	 * Overloaded < operator, compares if this container's error/numContributors < argument's error/numContributors
	 *
	 * @param container ErrorContainer to compare to
	 * @return true if container's error/numContributors < argument's error/numContributors
	 */
	bool operator< (const ErrorContainer<LT,IT>& container) const
	{
		return getContribPerPoint() < container.getContribPerPoint();
	}


	/**
	 * Returns the average error per contributor =  error/numContributors,
	 * or 0 if there are no contributors.
	 *
	 * @return average error per contributor =  error/numContributors
	 */
	RefinementFunctor::value_type getContribPerPoint() const
	{
		if(contributionCounter==0){
			return 0.0;
		}
		return error/static_cast<double>(contributionCounter);

		//		return error;
	}


	/**
	 * Returns if the grid point / subspace complies to the admissibility criterion of the
	 * algorithm
	 *
	 *@return true if the grid point / subspace complies to the admissibility criterion of the
	 *algorithm
	 */
	bool isAdmissible() const
	{
		return admissible;
	}

	/**
	 *Set the admissibility criterion for the grid point / subspace
	 *
	 *@param admissible set to true or false according
	 *to the admissibility criterion of the grid point / subspace
	 */
	void setAdmissible(bool admissible)
	{
		this->admissible = admissible;
	}

	/**
	 * Sets the errorIndicator and the contribution counter to 0.
	 */
	void resetError()
	{
		error = 0;
		contributionCounter = 0;
	}

	/**
	 *Outputs a string representation of the object with all its members.
	 *
	 *@return the object's string representation
	 */
	std::string toString()
	{

		std::ostringstream tmp;
		tmp << HashGridIndex<LT,IT>::toString()  << " with error : " << error << ", with: " << contributionCounter << " contributors ,mean: " << getContribPerPoint() << ". Refining this Object is admissible: " << admissible <<  "\n";
		return  tmp.str();
	}

private:

	// error indicator value
	RefinementFunctor::value_type error;
	// amount of contributors to the error indicator,
	//needed for a subspace indicator based on the mean
	size_t contributionCounter;
	// a flag to save if the grid point/subspace is admissible.
	bool admissible;

};

/**
 * class needed for hashing the underlying HashGridIndex class.
 */
template<class LT, class IT>
struct hash<ErrorContainer<LT, IT>* > {
	size_t operator()(ErrorContainer<LT, IT>* index) const {
		return index->hash();
	}
};

/**
 * class needed for the underlying HashGridIndex class to check two instances for equality.
 */
template<class LT, class IT>
struct eqIndex<ErrorContainer<LT, IT>* > {
	size_t operator()(ErrorContainer<LT, IT>* s1, ErrorContainer<LT, IT>* s2) const {
		return s1->equals(*s2);
	}
};

} /* namespace base */
} /* namespace SGPP */
#endif /* ErrorContainer_HPP_ */