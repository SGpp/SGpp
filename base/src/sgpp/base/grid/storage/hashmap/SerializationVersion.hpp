// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SERIALIZATIONVERSION_HPP
#define SERIALIZATIONVERSION_HPP

/**
 * This specifies the available serialization versions
 *
 * Version 1: classic verions without leaf proeperty
 * Version 2: every gridpoint is extended by one boolean that specifies if it's a leaf
 * Version 3: added support for the grid's bounding box
 * Version 4: needed for import of the Bonner's Sparse Grid Definition files; same as Ver 3
 *        but without leaf property, NOT FOR EXPORT
 * Version 5: differentiate BoundingBox and Stretching, added support for stretching.
 * Version 6: added PointDistribution to HashGridIndex
 *            ("Normal" with x = i*2^(-l) and "ClenshawCurtis")
 * Version 7: PointDistribution changed from enum to enum class
 * Version 8: Add custom boundaryLevel (>= 1) for LinearBoundaryGrid etc.
 */
#define SERIALIZATION_VERSION 8

#endif /* SERIALIZATIONVERSION_HPP */
