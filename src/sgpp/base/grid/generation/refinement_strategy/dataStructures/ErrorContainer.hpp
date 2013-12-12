/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef ErrorContainer_HPP_
#define ErrorContainer_HPP_

#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/storage/hashmap/HashGridIndex.hpp"
#include "base/grid/storage/hashmap/HashGridStorage.hpp"

namespace sg {
namespace base {

template<class LT,class IT>
class ErrorContainer: public HashGridIndex<LT,IT> {
public:

	ErrorContainer(RefinementFunctor::value_type error, size_t contributionCounter, HashGridIndex<LT,IT>& index):HashGridIndex<LT,IT>(index)
{
		this->error = error;
		this->contributionCounter = contributionCounter;
		admissible = true;
}

	ErrorContainer(RefinementFunctor::value_type error, HashGridIndex<LT,IT>& index):HashGridIndex<LT,IT>(index)
	{
		this->error = error;
		this->contributionCounter = 1;
		admissible = true;
	}

	ErrorContainer(HashGridIndex<LT,IT>& index):HashGridIndex<LT,IT>()
	{
		error = 0;
		contributionCounter = 0;
		admissible = true;
	}

	ErrorContainer(ErrorContainer<LT,IT>& container):HashGridIndex<LT,IT>(container)
	{
		error = container.error;
		contributionCounter = container.contributionCounter;
		admissible = container.admissible;
	}

	ErrorContainer(const ErrorContainer<LT,IT>& container):HashGridIndex<LT,IT>(container)
	{
		error = container.error;
		contributionCounter = container.contributionCounter;
		admissible = container.admissible;
	}

	ErrorContainer(ErrorContainer<LT,IT>* container):HashGridIndex<LT,IT>(container)
	{
		error = container->error;
		contributionCounter = container->contributionCounter;
		admissible = container->admissible;
	}

	ErrorContainer():HashGridIndex<LT,IT>()
	{
		error = 0;
		contributionCounter = 0;
		admissible = true;
	}

	RefinementFunctor::value_type getError()
	{
		return error;
	}

	size_t getContributionCounter()
	{
		return contributionCounter;
	}

	void operator= (RefinementFunctor::value_type newError)
	{
		error = newError;
		contributionCounter = 1;
	}


	ErrorContainer<LT,IT>& operator= (ErrorContainer<LT,IT>& container)
	{
		HashGridIndex<LT,IT>::operator =(container);
		error = container.error;
		contributionCounter = container.contributionCounter;
		admissible = container.admissible;

		return *this;
	}

	ErrorContainer<LT,IT> operator+ (RefinementFunctor::value_type newError)
	{
		ErrorContainer result = *this;
		result.error += newError;
		result.contributionCounter++;

		return result;
	}

	ErrorContainer<LT,IT> operator+ (ErrorContainer<LT,IT>& container)
	{
		ErrorContainer result = *this;
		result.error += container.error;
		result.contributionCounter += container.contributionCounter;

		return result;
	}

	ErrorContainer<LT,IT> operator+= (RefinementFunctor::value_type newError)
	{
		this->error += newError;
		this->contributionCounter ++;
		return *this;
	}

	ErrorContainer<LT,IT>& operator+= (ErrorContainer<LT,IT>& container)
	{
		this->error += container.error;
		this->contributionCounter += container.contributionCounter;
		return *this;
	}

	ErrorContainer<LT,IT> operator- (RefinementFunctor::value_type newError)
	{
		ErrorContainer result = *this;
		result.error -= newError;
		result.contributionCounter++;

		return result;
	}


	ErrorContainer<LT,IT> operator- (ErrorContainer<LT,IT>& container)
	{
		ErrorContainer result = *this;
		result.error -= container.error;
		result.contributionCounter += container.contributionCounter;

		return result;
	}

	bool operator== (RefinementFunctor::value_type otherValue)
	{
		return getContribPerPoint() == otherValue;
	}

	bool operator== (ErrorContainer<LT,IT>& container)
	{
		return getContribPerPoint() == container.getContribPerPoint();
	}

	bool operator> (RefinementFunctor::value_type otherValue)
	{
		return getContribPerPoint() > otherValue;
	}

	bool operator> (ErrorContainer<LT,IT>& container)
	{
		return getContribPerPoint() > container.getContribPerPoint();
	}

	bool operator> (const ErrorContainer<LT,IT>& container) const
	{
		return getContribPerPoint() > container.getContribPerPoint();
	}

	bool operator< (RefinementFunctor::value_type otherValue)
	{
		return getContribPerPoint() < otherValue;
	}

	bool operator< (ErrorContainer<LT,IT>& container)
	{
		return getContribPerPoint() < container.getContribPerPoint();
	}

	bool operator< (const ErrorContainer<LT,IT>& container) const
	{
		return getContribPerPoint() < container.getContribPerPoint();
	}

	RefinementFunctor::value_type getContribPerPoint()
	{
		if(contributionCounter==0){
			return 0.0;
		}
		return error/static_cast<double>(contributionCounter);
	}

	RefinementFunctor::value_type getContribPerPoint() const
	{
		if(contributionCounter==0){
			return 0.0;
		}
		return error/static_cast<double>(contributionCounter);
	}


	bool isAdmissible() const
	{
		return admissible;
	}

	void setAdmissible(bool admissible)
	{
		this->admissible = admissible;
	}

	std::string toString()
	{

		std::ostringstream tmp;
		tmp << HashGridIndex<LT,IT>::toString()  << " with error : " << error << ", with: " << contributionCounter << " contributors ,mean: " << getContribPerPoint() << ". Refining this Object is admissible: " << admissible <<  "\n";
		return  tmp.str();
	}

private:

	RefinementFunctor::value_type error;
	size_t contributionCounter;
	bool admissible;


};


template<class LT, class IT>
    struct hash<ErrorContainer<LT, IT>* > {
      size_t operator()(ErrorContainer<LT, IT>* index) const {
        return index->hash();
      }
    };

    template<class LT, class IT>
    struct eqIndex<ErrorContainer<LT, IT>* > {
      size_t operator()(ErrorContainer<LT, IT>* s1, ErrorContainer<LT, IT>* s2) const {
        return s1->equals(*s2);
      }
    };


} /* namespace base */
} /* namespace sg */
#endif /* ErrorContainer_HPP_ */
