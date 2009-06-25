/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2007 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef DATAVECTORDEFINITION
#define DATAVECTORDEFINITION

/**
 * This struct is needed for exporting a DataVector
 * to another address space, so it contains all 
 * information that is needed to reconstruct a 
 * DataVector object
 *
 * The spave required by a DataVector object is:
 * (size+unused)*dim*sizeof(double)
 */
struct DataVectorDefinition
{
	/// Pointer to the data of the DataVector
	double* pointerToData;
	/// Number of Dimensions
	int dim;
	/// Number of Elements per Dim
	int size;
	/// Number of unused slots;
	int unused;
};

#endif /* DATAVECTORDEFINITION */
