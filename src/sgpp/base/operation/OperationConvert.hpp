/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
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

#ifndef OPERATIONCONVERT_HPP
#define OPERATIONCONVERT_HPP

#include "base/datatypes/DataVector.hpp"

#ifdef _WIN32
#pragma warning(disable: 4267)
#endif

namespace sg
{
namespace base
{

/**
 * Operation that cpnverts a given basis into the normal, linear hat basis and vice versa
 *
 */
class OperationConvert
{
public:
	/**
	 * Constructor
	 */
	OperationConvert() {}

	/**
	 * Destructor
	 */
	virtual ~OperationConvert() {}

	/**
	 * Convert given basis into linear hat basis.
	 */
	virtual void doConvertToLinear(DataVector& alpha) = 0;


	/**
	 * Convert from a linear coefficient vector into given basis.
	 */
	virtual void doConvertFromLinear(DataVector& alpha) = 0;
};

}
}

#endif /* OPERATIONCONVERT_HPP */
