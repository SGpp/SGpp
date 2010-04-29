/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef GENERATION_EXCEPTION_HPP
#define GENERATION_EXCEPTION_HPP

#include <exception>

namespace sg
{

/**
 * Exception that is thrown in case of a grid generation failure
 *
 * @version $HEAD$
 */
class generation_exception : public std::exception
{
public:
	/**
	 * Constructor
	 *
	 * @param msg the exception message
	 */
	generation_exception(const char* msg) throw() : msg(msg)
	{
	}

	/**
	 * Standared Constructor
	 */
	generation_exception() throw() : msg(NULL) { }

	/**
	 * Destructor
	 */
    virtual ~generation_exception() throw() { }

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
			return "generation_exception: failure generating grid";
		}
	}

protected:
	/// the exception message
	const char* msg;
};

}

#endif /* GENERATION_EXCEPTION_HPP */
