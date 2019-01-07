// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/*
 * AnisotropicFullGrid.cpp
 *
 *  Created on: Jan 5, 2019
 *      Author: Nico Roesel
 */

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/AnisotropicFullGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>
#include <vector>


namespace sgpp {
namespace base {

AnisotropicFullGrid::AnisotropicFullGrid(std::istream& istr) :
		Grid(istr),
		generator(storage){
}

AnisotropicFullGrid::AnisotropicFullGrid(size_t dim, std::vector<size_t> v) :
	Grid(dim),
	generator(storage){
}

AnisotropicFullGrid::~AnisotropicFullGrid() {
}

sgpp::base::GridType AnisotropicFullGrid::getType() {
	return sgpp::base::GridType::AnisotropicFullGrid;
}

SBasis& AnisotropicFullGrid::getBasis() {
	static SLinearBase basis;
	return basis;
}

GridGenerator& AnisotropicFullGrid::getGenerator(){
	return generator;
}

Grid* AnisotropicFullGrid::unserialize(std::istream& istr) {
	return new AnisotropicFullGrid(istr);
}

}
}


