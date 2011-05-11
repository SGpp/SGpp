/*
 * TruncatedTrapezoidGridGenerator.hpp
 *
 *  Created on: Aug 4, 2010
 *      Author: Aliz Nagy
 */

#ifndef TRUNCATEDTRAPEZOIDGRIDGENERATOR_HPP_
#define TRUNCATEDTRAPEZOIDGRIDGENERATOR_HPP_
#include "grid/GridStorage.hpp"
#include "grid/generation/GridGenerator.hpp"

namespace sg
{
namespace base
{

/**
 * This class provides the interface for the grid generation
 * for grids with boundaries, pentagon cut through sub space scheme
 */
class TruncatedTrapezoidGridGenerator : public GridGenerator
{
public:
	/**
	 * Constructor
	 *
	 * @param storage template type that holds the grid points
	 */
	TruncatedTrapezoidGridGenerator(GridStorage* storage);

	/**
	 * Destructor
	 */
	virtual ~TruncatedTrapezoidGridGenerator();
	/**
	 * Creates a regular trapezoid boundary grid with given level and l_user=1
	 * Is the same as the regular trapezoid grid
	 * */
	virtual void regular(size_t level);
	/**
	 * Creates a super trapezoid boundary grid with given level and l_user
	 * @param level the maximum level of the grid
	 * @param l_user the number of fullgrids cut off from the boundaries.
	 * */
	virtual void truncated(size_t level,size_t l_user);
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
}

#endif /* TRUNCATEDTRAPEZOIDGRIDGENERATOR_HPP_ */
