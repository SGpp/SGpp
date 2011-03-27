/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TOOL_EXCEPTION_HPP
#define TOOL_EXCEPTION_HPP

#include <exception>

namespace sg
{

/**
 * Exception that is thrown in case of a tool operation failure
 *
 * @version $HEAD$
 */
class tool_exception : public std::exception
{
public:
	/**
	 * Constructor
	 *
	 * @param msg the exception message
	 */
	tool_exception(const char* msg) throw() : msg(msg)
	{
	}

	/**
	 * Standard Constructor
	 */
	tool_exception() throw() : msg(NULL) { }

	/**
	 * Destructor
	 */
    virtual ~tool_exception() throw() { }

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
			return "tool_exception: general failure";
		}
	}

protected:
	/// the exception message
	const char* msg;

};

}

#endif /* TOOL_EXCEPTION_HPP */
