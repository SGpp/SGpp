/*
 * SquareRootGridGenerator.hpp
 *
 *  Created on: Aug 4, 2010
 *      Author: Aliz Nagy
 */

#ifndef SQUAREROOTGRIDGENERATOR_HPP_
#define SQUAREROOTGRIDGENERATOR_HPP_
#include "grid/GridStorage.hpp"
#include "grid/generation/GridGenerator.hpp"

namespace sg
{

/**
 * This class provides the interface for the grid generation
 * for grids with boundaries, pentagon cut through sub space scheme
 */
class SquareRootGridGenerator : public GridGenerator
{
public:
	/**
	 * Constructor
	 *
	 * @param storage template type that holds the grid points
	 */
	SquareRootGridGenerator(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~SquareRootGridGenerator();

	virtual void regular(size_t level);
	virtual void refine(RefinementFunctor* func){};
	virtual size_t getNumberOfRefinablePoints(){return 0;};

	virtual void coarsen(CoarseningFunctor* func, DataVector* alpha){};
	virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly){};
	virtual size_t getNumberOfRemoveablePoints(){return 0;};

	virtual void refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel){};
	virtual size_t getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel){return 0;};

protected:
	/// Pointer to the grid's storage object
	GridStorage* storage;
};

}

#endif /* SQUAREROOTGRIDGENERATOR_HPP_ */
