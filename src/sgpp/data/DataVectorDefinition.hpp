/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DATAVECTORDEFINITION_HPP
#define DATAVECTORDEFINITION_HPP

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

#endif /* DATAVECTORDEFINITION_HPP */
