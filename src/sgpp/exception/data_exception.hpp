/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (Dirk.Pflueger@in.tum.de)

#ifndef DATA_EXCEPTION_HPP
#define DATA_EXCEPTION_HPP

#include <exception>

namespace sg
{

/**
 * Exception that is thrown in case of a data failure (conversion, creation, ...).
 *
 * @version $HEAD$
 */
class data_exception : public std::exception
{
public:
	/**
	 * Create a new exception (constructor) with some message.
	 *
	 * @param msg The exception message
	 */
	data_exception(const char* msg) throw() : msg(msg)
	{
	}

	/**
	 * Create default exception (constructor).
	 */
	data_exception() throw() : msg(NULL) { }

	/**
	 * Destructor
	 */
    virtual ~data_exception() throw() { }

    /**
     * Return message of exception object.
     *
     * @return Returns the message specified in the constructor, otherwise a general text
     */
	virtual const char* what() const throw()
	{
		if(msg)
		{
			return msg;
		}
		else
		{
			return "data_exception: general failure";
		}
	}

protected:
	/// the exception message
	const char* msg;

};

}

#endif /* DATA_EXCEPTION_HPP */
