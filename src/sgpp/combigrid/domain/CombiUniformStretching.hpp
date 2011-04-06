/*
 * CombiUniformStretching.hpp
 *
 *  Created on: Apr 4, 2011
 *      Author: benk
 */

#ifndef COMBIUNIFORMSTRETCHING_HPP_
#define COMBIUNIFORMSTRETCHING_HPP_

#include "combigrid/domain/AbstractStretchingMaker.hpp"

namespace combigrid {

/** uniform stretching used only for testing purposes */
class UniformStretching : public AbstractStretchingMaker{
public:

	UniformStretching():AbstractStretchingMaker() {;}

	virtual ~UniformStretching() {;}

	void get1DStretching(
			int level , double min, double max,
			std::vector<double>& stretching) const;

};

}

#endif /* COMBIUNIFORMSTRETCHING_HPP_ */
