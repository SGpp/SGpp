/*
 * StretchingFactory.cpp
 *
 *  Created on: Jan 19, 2015
 *      Author: petz
 */

#include <sgpp/combigrid/domain/StretchingFactory.hpp>

namespace combigrid {


AbstractStretchingMaker* createStretchingMaker(combigrid::Stretching str) {

		AbstractStretchingMaker* maker = NULL;

		switch (str) {
		case ATAN:
			maker = new AtanSpecialStretching();
			break;
		case TAN:
			maker = new TanStretching();
			break;
		case CHEBYSHEV:
			maker = new CombiChebyshevStretching();
			break;
		case EQUIDISTANT:
			maker = new CombiEquidistantStretching();
			break;
		case LEGENDRE:
			maker = new CombiLegendreStretching();
			break;
		case BASU:
			maker = new CombiBasuStretching();
			break;
		case UNKNOWN:
		default:
			COMBIGRID_OUT_WRN("Input stretching type unknown! Returning a NULL stretching maker!",__FILE__,__LINE__)
			break;
		}

		return maker;

	}

}
