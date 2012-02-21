/*
 * AdaptiveSerialCombigridVariableCoefficients.hpp
 *
 *  Created on: 21.02.2012
 *      Author: ckow
 */

#ifndef ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_
#define ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_

#include "SerialCombiGrid.hpp"

namespace combigrid {

class AdaptiveSerialCombigridVariableCoefficients: public combigrid::SerialCombiGrid {
public:
	AdaptiveSerialCombigridVariableCoefficients(
			const CombiSchemeBasis * combischeme,
			const std::vector<bool>& hasBoundaryPts) :
			SerialCombiGrid(combischeme, hasBoundaryPts) {
		;
	}
	virtual ~AdaptiveSerialCombigridVariableCoefficients();

	void changeCoefficient();
};

} /* namespace combigrid */
#endif /* ADAPTIVESERIALCOMBIGRIDVARIABLECOEFFICIENTS_HPP_ */
