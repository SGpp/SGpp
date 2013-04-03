/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Richard Röttger

#ifndef PREWAVELETGRIDGENERATOR_HPP
#define PREWAVELETGRIDGENERATOR_HPP

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/GridGenerator.hpp"

namespace sg
{
namespace base
{

/**
 * GridGenerator for prewavelet grids without boundaries
 */
class PrewaveletGridGenerator: public GridGenerator
{
protected:
	/// pointer to the storage object
	GridStorage* storage;
	GridStorage* shadowstorage;
	typedef GridStorage::index_type index_type;
	typedef index_type::index_type index_t;
	typedef index_type::level_type level_t;

private:
	void addNeighbours(index_type& index, size_t current_dim,
			level_t target_level, GridStorage::grid_iterator& iter, GridStorage::grid_iterator& shadowIter);
	void insertParents(GridStorage::grid_iterator& iter, GridStorage::grid_iterator& shadowIter);
	void consolidateShadow();

public:
	/**
	 * Constructor
	 *
	 * @param storage the grid's storage object
	 */
	PrewaveletGridGenerator(GridStorage* storage, GridStorage* shadowstorage);

	/**
	 * Destructor
	 */
	virtual ~PrewaveletGridGenerator();

	virtual void regular(size_t level);
	virtual void full(size_t level);
	virtual void refine(RefinementFunctor* func);
	virtual size_t getNumberOfRefinablePoints();

	virtual void coarsen(CoarseningFunctor* func, DataVector* alpha);
	virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha, size_t numFirstOnly);
	virtual size_t getNumberOfRemovablePoints();

	virtual void refineMaxLevel(RefinementFunctor* func, unsigned int maxLevel);
	virtual size_t getNumberOfRefinablePointsToMaxLevel(unsigned int maxLevel);

};

}
}

#endif /* PREWAVELETGRIDGENERATOR_HPP */
