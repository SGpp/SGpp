/*****************************************************************************/
/* This file is part of sg++, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sg++ is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sg++ is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sg++; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef OPERATIONEVALMODWAVELET_HPP
#define OPERATIINEVALMODWAVELET_HPP

#include "operation/OperationEval.hpp"
#include "grid/GridStorage.hpp"

namespace sg
{

class OperationEvalModWavelet : public OperationEval
{
public:
	OperationEvalModWavelet(GridStorage* storage) : storage(storage) {}
	virtual ~OperationEvalModWavelet() {}

	virtual double eval(DataVector& alpha, std::vector<double>& point);
	virtual double test(DataVector& alpha, DataVector& data, DataVector& classes);

protected:
	GridStorage* storage;

};

}

#endif /* OPERATIINEVALMODWAVELET_HPP */
