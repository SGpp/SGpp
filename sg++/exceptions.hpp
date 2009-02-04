/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2008 Joerg Blank (blankj@in.tum.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EXCEPTIONS_HPP_
#define EXCEPTIONS_HPP_

#include <exception>


namespace sg
{

class generation_exception : public std::exception
{
public:
	generation_exception(const char* msg) throw() : msg(msg)
	{
	}

	generation_exception() throw() : msg(NULL) { }

    virtual ~generation_exception() throw() { }

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
	const char* msg;

};

class factory_exception : public std::exception
{
public:
	factory_exception(const char* msg) throw() : msg(msg)
	{
	}

	factory_exception() throw() : msg(NULL) { }

    virtual ~factory_exception() throw() { }

	virtual const char* what() const throw()
	{
		if(msg)
		{
			return msg;
		}
		else
		{
			return "factory_exception: general failure";
		}
	}
protected:
	const char* msg;

};

class operation_exception : public std::exception
{
public:
	operation_exception(const char* msg) throw() : msg(msg)
	{
	}

	operation_exception() throw() : msg(NULL) { }

    virtual ~operation_exception() throw() { }

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
	const char* msg;

};

}

#endif /*EXCEPTIONS_HPP_*/
