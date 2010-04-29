/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef ALGORITHM_EXCEPTION_HPP
#define ALGORITHM_EXCEPTION_HPP

#include <exception>

namespace sg
{

/**
 * Exception that is thrown in case of an application failure
 *
 * @version $HEAD$
 */
class algorithm_exception : public std::exception
{
public:
	/**
	 * Constructor
	 *
	 * @param msg the exception message
	 */
	algorithm_exception(const char* msg) throw() : msg(msg)
	{
	}

	/**
	 * Standard Constructor
	 */
	algorithm_exception() throw() : msg(NULL) { }

	/**
	 * Destructor
	 */
    virtual ~algorithm_exception() throw() { }

    /**
     * throw method that have to be implemented
     *
     * @return returns the message specified in the constructor otherwise a general text
     */
	virtual const char* what() const throw()
	{
		if(msg)
		{
			return msg;
		}
		else
		{
			return "algorithm_exception: general failure";
		}
	}

protected:
	/// the exception message
	const char* msg;

};

}

#endif /* ALGORITHM_EXCEPTION_HPP */
