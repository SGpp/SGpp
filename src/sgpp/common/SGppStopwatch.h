/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef SGPPSTOPWATCH_H
#define SGPPSTOPWATCH_H

#ifdef WINDOWS
#include <windows.h>
#endif
#ifndef WINDOWS
#include <sys/time.h>
#endif

namespace sg
{

/**
 *	OS-independent (per Preprocessor) version of a stopwatch
 *
 *	Part of SGpp, so you can easily calculate the needed time of SGpp computations with a high precision
 */
class SGppStopwatch
{
private:
#ifdef WINDOWS
	LARGE_INTEGER ticksPerSecond;
	LARGE_INTEGER begin;
	LARGE_INTEGER end;
#endif
#ifndef WINDOWS
	timeval begin;
	timeval end;
#endif

public:
	/**
	 *	Constructor
	 *
	 *	resets the Stopwatch
	 */
	SGppStopwatch();

	/**
	 *	Destructor
	 */
	~SGppStopwatch();

	/**
	 *	starts the time measuring
	 */
	void start();

	/**
	 *	stops time measuring
	 *
	 *	\return measured time in seconds
	 */
	double stop();
};

}

#endif	/* SGPPSTOPWATCH_H */
