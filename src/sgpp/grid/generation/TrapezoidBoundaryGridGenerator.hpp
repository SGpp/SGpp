/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TRAPEZOIDBOUNDARYGRIDGENERATOR_HPP
#define TRAPEZOIDBOUNDARYGRIDGENERATOR_HPP

#include "grid/GridStorage.hpp"
#include "grid/generation/GridGenerator.hpp"

namespace sg
{

/**
 * This class provides the interface for the grid generation
 * for grids with boundaries, pentagon cut through sub space scheme
 */
class TrapezoidBoundaryGridGenerator : public GridGenerator
{
public:
	/**
	 * Constructor
	 *
	 * @param storage template type that holds the grid points
	 */
	TrapezoidBoundaryGridGenerator(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~TrapezoidBoundaryGridGenerator();

	virtual void regular(size_t level);
	virtual void refine(RefinementFunctor* func);
	virtual int getNumberOfRefinablePoints();

protected:
	/// Pointer to the grid's storage object
	GridStorage* storage;
};

}

#endif /* TRAPEZOIDBOUNDARYGRIDGENERATOR_HPP */
