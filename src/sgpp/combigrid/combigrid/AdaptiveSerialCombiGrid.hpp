/*
 * AdaptiveSerialCombiGrid.hpp
 *
 *  Created on: Jun 28, 2011
 *      Author: kowitz_local
 */

#ifndef ADAPTIVESERIALCOMBIGRID_HPP_
#define ADAPTIVESERIALCOMBIGRID_HPP_

#include "combigrid/combigrid/SerialCombiGrid.hpp"


namespace combigrid {
/**
 * Class allowing the dimension adaptive change of the combischeme
 */
class AdaptiveSerialCombiGrid: public SerialCombiGrid {
public:
	/**
	 * Constructor similar requiring the initial combischeme
	 */
	AdaptiveSerialCombiGrid(const CombiSchemeBasis * combischeme,
			const std::vector<bool>& hasBoundaryPts) :
		SerialCombiGrid(combischeme, hasBoundaryPts) {
		;
	}

//	~AdaptiveSerialCombiGrid();

	/** Method to add an extra fullgrid to the combination.
	 * The coefficients are set accordingly ( zero coefficients are allowed).
	 * returns the indices of fullgrids, which have to be filled with new data
	 */
	std::vector<int> addToCombiScheme(std::vector<int> level);
};

}

#endif /* ADAPTIVESERIALCOMBIGRID_HPP_ */
