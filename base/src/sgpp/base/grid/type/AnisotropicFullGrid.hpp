// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/*
 * AnisotropicFullGrid.hpp
 *
 *  Created on: Jan 5, 2019
 *      Author: Nico Roesel
 */

#ifndef ANISOTROPICFULLGRID_HPP
#define ANISOTROPICFULLGRID_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {
/**
* Anisotropic full grid. Used for example in the sparse grid combination technique.
*/

class AnisotropicFullGrid : public Grid {
protected:
	///grid generator
	StandardGridGenerator generator;
	explicit AnisotropicFullGrid(std::istream& istr);

public:
	/**
	 * Constructor Anisotropic Full Grid
	 *
	 * @param dim the dimension of the grid
	 * @param v the level vector of the dimensions,
	 * its size must equal dim
	 */
	explicit AnisotropicFullGrid(size_t dim, std::vector<size_t> v);


	/**
	 * Destructor
	 */
	~AnisotropicFullGrid() override;

	sgpp::base::GridType getType() override;

	SBasis& getBasis() override;

	GridGenerator& getGenerator() override;

	static Grid* unserialize(std::istream& istr);
};
}}



#endif /* ANISOTROPICFULLGRID_HPP */
