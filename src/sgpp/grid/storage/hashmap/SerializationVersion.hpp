/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SERIALIZATIONVERSION_HPP
#define SERIALIZATIONVERSION_HPP

/**
 * This specifies the available serialization versions
 *
 * Version 1: classic verions without leaf proeperty
 * Version 2: every gridpoint is extended by one boolean that specifies if it's a leaf
 * Version 3: added support for the grid's bounding box
 * Version 4: needed for import of the Bonner's Sparse Grid Definition files; same as Ver 3
 * 			  but without leaf property, NOT FOR EXPORT
 */
#define SERIALIZATION_VERSION 3

#endif /* SERIALIZATIONVERSION_HPP */
