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

#include "tools/common/SGppStopwatch.hpp"

namespace sg
{

SGppStopwatch::SGppStopwatch()
{
#ifdef WINDOWS
	QueryPerformanceFrequency(&ticksPerSecond);
#endif
#ifndef WINDOWS

#endif
}

SGppStopwatch::~SGppStopwatch()
{
}

void SGppStopwatch::start()
{
#ifdef WINDOWS
	QueryPerformanceCounter(&begin);
#endif
#ifndef WINDOWS
	gettimeofday(&begin,(struct timezone *)0);
#endif
}

double SGppStopwatch::stop()
{
#ifdef WINDOWS
	QueryPerformanceCounter(&end);

	double ret, ticksps;

	end.QuadPart -= begin.QuadPart;
	ret = (double)(end.QuadPart);
	ticksps = (double)(ticksPerSecond.QuadPart);
	ret /= ticksps;

	return ret;
#endif
#ifndef WINDOWS
	gettimeofday(&end,(struct timezone *)0);
	double seconds, useconds;
	double ret, tmp;

	if (end.tv_usec >= begin.tv_usec)
	{
		seconds = (double)end.tv_sec - (double)begin.tv_sec;
		useconds = (double)end.tv_usec - (double)begin.tv_usec;
	}
	else
	{
		seconds = (double)end.tv_sec - (double)begin.tv_sec;
		seconds -= 1;					// Correction
		useconds = (double)end.tv_usec - (double)begin.tv_usec;
		useconds += 1000000;			// Correction
	}

	// get time in seconds
	tmp = (double)useconds;
	ret = (double)seconds;
	tmp /= 1000000;
	ret += tmp;

	return ret;
#endif
}

}
