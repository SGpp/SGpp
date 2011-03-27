/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATION_EXCEPTION_HPP
#define OPERATION_EXCEPTION_HPP

#include <exception>

namespace sg
{

/**
 * Exception that is thrown in case of a grid operation failure
 *
 * @version $HEAD$
 */
class operation_exception : public std::exception
{
public:
	/**
	 * Constructor
	 *
	 * @param msg the exception message
	 */
	operation_exception(const char* msg) throw() : msg(msg)
	{
	}

	/**
	 * Standard Constructor
	 */
	operation_exception() throw() : msg(NULL) { }

	/**
	 * Destructor
	 */
    virtual ~operation_exception() throw() { }

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
			return "operation_exception: general failure";
		}
	}

protected:
	/// the exception message
	const char* msg;

};

}

#endif /* OPERATION_EXCEPTION_HPP */
