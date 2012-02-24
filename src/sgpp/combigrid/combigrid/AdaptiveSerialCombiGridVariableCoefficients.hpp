/*
 * AdaptiveSerialCombiGridVariableCoefficients.hpp
 *
 *  Created on: 21.02.2012
 *      Author: ckow
 */

#ifndef ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_
#define ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_

#include "combigrid/combigrid/SerialCombiGrid.hpp"
#include "combigrid/utils/combigrid_ultils.hpp"

namespace combigrid {

class AdaptiveSerialCombiGridVariableCoefficients: public combigrid::SerialCombiGrid {
public:
	AdaptiveSerialCombiGridVariableCoefficients(
			const CombiSchemeBasis * combischeme,
			const std::vector<bool>& hasBoundaryPts) :
			SerialCombiGrid(combischeme, hasBoundaryPts) {
		;
	}
//	virtual ~AdaptiveSerialCombiGridVariableCoefficients();

	void changeCoefficients(std::vector<double> newCoef);
	void changeCoefficients(int i,double newCoef);
};

} /* namespace combigrid */
#endif /* ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_ */
