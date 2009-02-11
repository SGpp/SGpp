/*****************************************************************************/
/* This file is part of pysgpp, a program package making use of spatially    */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* pysgpp is free software; you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* pysgpp is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with pysgpp; if not, write to the Free Software                     */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

namespace sg
{

class RefinementFunctor
{
public:
	typedef double value_type;

	virtual double operator()(GridStorage* storage, size_t seq) = 0;
	virtual double start() = 0;	
};

class GridGenerator
{
public:
	virtual void regular(size_t level) = 0;
	virtual void refine(RefinementFunctor* func) = 0;
};

class OperationB
{
public:
	virtual void mult(DataVector& alpha, DataVector& data, DataVector& result) = 0;
	virtual void multTranspose(DataVector& alpha, DataVector& data, DataVector& result) = 0;
};

class OperationMatrix
{
public:
	virtual void mult(DataVector& alpha, DataVector& result) = 0;
};

class OperationEval
{
public:
	virtual double eval(DataVector& alpha, DataVector& point) = 0;
	virtual double test(DataVector& alpha, DataVector& data, DataVector& classes) = 0;
};

class OperationHierarchisation
{
public:
	virtual void doHierarchisation(DataVector& node_values) = 0;
	virtual void doDehierarchisation(DataVector& alpha) = 0;
};

}
