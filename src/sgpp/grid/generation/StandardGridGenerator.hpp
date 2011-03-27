/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef STANDARDGRIDGENERATOR_HPP
#define STANDARDGRIDGENERATOR_HPP

#include "grid/GridStorage.hpp"
#include "grid/generation/GridGenerator.hpp"

namespace sg
{

/**
 * GridGenerator for standard grids without boundaries
 */
class StandardGridGenerator : public GridGenerator
{
public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's storage object
	 */
	StandardGridGenerator(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~StandardGridGenerator();

	virtual void regular(size_t level);
	virtual void refine(RefinementFunctor* func);
	virtual int getNumberOfRefinablePoints();

protected:
	/// pointer to the storage object
	GridStorage* storage;
};

}

#endif /* STANDARDGRIDGEMERATOR_HPP */
