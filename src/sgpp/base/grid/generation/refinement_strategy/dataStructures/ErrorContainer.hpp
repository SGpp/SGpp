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

namespace sg {
namespace base {

class ErrorContainer {
public:

	ErrorContainer(RefinementFunctor::value_type error, size_t contributionCounter);

	ErrorContainer(RefinementFunctor::value_type error);

	ErrorContainer();

	RefinementFunctor::value_type getError();

	size_t getContributionCounter();

	void operator= (RefinementFunctor::value_type newError);

	ErrorContainer operator+ (RefinementFunctor::value_type newError);

	ErrorContainer operator+ (ErrorContainer& container);

	ErrorContainer operator+= (RefinementFunctor::value_type newError);

	ErrorContainer& operator+= (ErrorContainer& container);

	ErrorContainer operator- (RefinementFunctor::value_type newError);

	ErrorContainer operator- (ErrorContainer& container);

	bool operator== (RefinementFunctor::value_type otherValue);

	bool operator== (ErrorContainer& container);

	bool operator> (RefinementFunctor::value_type otherValue);

	bool operator> (ErrorContainer& container);

	bool operator< (RefinementFunctor::value_type otherValue);

	bool operator< (ErrorContainer& container);

	RefinementFunctor::value_type getContribPerPoint();

	std::string toString();

private:

	RefinementFunctor::value_type error;
	size_t contributionCounter;


};

} /* namespace base */
} /* namespace sg */
#endif /* ErrorContainer_HPP_ */
