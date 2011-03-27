/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BOUNDARYGRIDGENERATOR_HPP
#define BOUNDARYGRIDGENERATOR_HPP

#include "grid/GridStorage.hpp"
#include "grid/generation/GridGenerator.hpp"

namespace sg
{

/**
 * This class provides the interface for the grid generation
 * for grids with boundaries, diagonal cut through sub space scheme
 */
class BoundaryGridGenerator : public GridGenerator
{
public:
	/**
	 * Constructor
	 *
	 * @param storage template type that holds the grid points
	 */
	BoundaryGridGenerator(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~BoundaryGridGenerator();

	virtual void regular(size_t level);
	virtual void refine(RefinementFunctor* func);
	virtual int getNumberOfRefinablePoints();

protected:
	/// Pointer to the grid's storage object
	GridStorage* storage;
};

}

#endif /* BOUNDARYGRIDGEMERATOR_HPP */
