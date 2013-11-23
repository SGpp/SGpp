/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#include "ErrorContainer.hpp"

namespace sg {
namespace base {

ErrorContainer::ErrorContainer(RefinementFunctor::value_type error,
		size_t contributionCounter) {
	this->error = error;
	this->contributionCounter = contributionCounter;
}

ErrorContainer::ErrorContainer(RefinementFunctor::value_type error) {
	this->error = error;
	this->contributionCounter = 1;
}

ErrorContainer::ErrorContainer() {
	error = 0;
	contributionCounter = 0;
}

RefinementFunctor::value_type ErrorContainer::getError() {
	return error;
}

size_t ErrorContainer::getContributionCounter() {
	return contributionCounter;
}

void ErrorContainer::operator =(
		RefinementFunctor::value_type newError) {
	error = newError;
	contributionCounter = 1;
}

ErrorContainer ErrorContainer::operator +(
		RefinementFunctor::value_type newError) {

	ErrorContainer result = *this;
	result.error += newError;
	result.contributionCounter++;

	return result;
}

ErrorContainer ErrorContainer::operator +(
		ErrorContainer& container) {
	ErrorContainer result = *this;
	result.error += container.error;
	result.contributionCounter += container.contributionCounter;

	return result;
}

ErrorContainer ErrorContainer::operator +=(
		RefinementFunctor::value_type newError) {
	this->error += newError;
	this->contributionCounter ++;
	return *this;
}

ErrorContainer& ErrorContainer::operator +=(
		ErrorContainer& container) {
	this->error += container.error;
	this->contributionCounter += container.contributionCounter;
	return *this;
}

ErrorContainer ErrorContainer::operator -(
		RefinementFunctor::value_type newError) {

	ErrorContainer result = *this;
	result.error -= newError;
	result.contributionCounter++;

	return result;
}

ErrorContainer ErrorContainer::operator -(
		ErrorContainer& container) {

	ErrorContainer result = *this;
	result.error -= container.error;
	result.contributionCounter += container.contributionCounter;

	return result;
}

bool ErrorContainer::operator ==(
		RefinementFunctor::value_type otherValue) {
	return error == otherValue ? true : false ;
}

bool ErrorContainer::operator ==(ErrorContainer& container) {
	return error == container.error ? true : false ;
}

bool ErrorContainer::operator >(
		RefinementFunctor::value_type otherValue) {
	return error > otherValue ? true : false ;
}

bool ErrorContainer::operator >(ErrorContainer& container) {
	return error > container.error ? true : false ;
}

bool ErrorContainer::operator <(
		RefinementFunctor::value_type otherValue) {
	return error < otherValue ? true : false ;
}

bool ErrorContainer::operator <(ErrorContainer& container) {
	return error < container.error ? true : false ;
}

RefinementFunctor::value_type ErrorContainer::getContribPerPoint() {
	return error/static_cast<double>(contributionCounter);
}

std::string ErrorContainer::toString() {
	std::ostringstream tmp;
	tmp << "error: " << error << ", iterations: " << contributionCounter << ",contrib: " << getContribPerPoint() <<  "\n";
	return  tmp.str();
}


} /* namespace base */
} /* namespace sg */

