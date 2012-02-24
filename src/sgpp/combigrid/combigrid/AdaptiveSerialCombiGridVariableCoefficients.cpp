/*
 * AdaptiveSerialCombiGridVariableCoefficients.cpp
 *
 *  Created on: 21.02.2012
 *      Author: ckow
 */

#include "combigrid/combigrid/AdaptiveSerialCombiGridVariableCoefficients.hpp"
/* namespace combigrid */

void combigrid::AdaptiveSerialCombiGridVariableCoefficients::changeCoefficients(
		std::vector<double> newCoef) {
//	COMBIGRID_OUT("Setting new coefficients in combischeme.");

	combischeme_->setCoef(newCoef);
	combikernel_->setCoef(newCoef);
}

void combigrid::AdaptiveSerialCombiGridVariableCoefficients::changeCoefficients(
		int i, double newCoef) {
//	COMBIGRID_OUT("Setting new coefficient in combischeme.");

	combischeme_->setCoef(i,newCoef);
	combikernel_->setCoef(i,newCoef);
}
