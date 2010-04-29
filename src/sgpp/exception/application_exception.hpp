/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef APPLICATION_EXCEPTION_HPP
#define APPLICATION_EXCEPTION_HPP

#include <exception>

namespace sg
{

/**
 * Exception that is thrown in case of an application failure
 *
 * @version $HEAD$
 */
class application_exception : public std::exception
{
public:
	/**
	 * Constructor
	 *
	 * @param msg the exception message
	 */
	application_exception(const char* msg) throw() : msg(msg)
	{
	}

	/**
	 * Standard Constructor
	 */
	application_exception() throw() : msg(NULL) { }

	/**
	 * Destructor
	 */
    virtual ~application_exception() throw() { }

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
			return "application_exception: general failure";
		}
	}

protected:
	/// the exception message
	const char* msg;

};

}

#endif /* APPLICATION_EXCEPTION_HPP */
