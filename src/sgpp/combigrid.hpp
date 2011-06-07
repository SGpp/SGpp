/*
 * combigrid.hpp
 *
 *  Created on: Mar 23, 2011
 *      Author: benk
 */

#ifndef COMBIGRID_HPP_
#define COMBIGRID_HPP_

// --------- Core SGpp includes -----------
#include "data/DataVector.hpp"
#include "grid/Grid.hpp"
#include "grid/GridStorage.hpp"

// ----------- Combi Grid Includes ----------
#include "combigrid/combigrid/AbstractCombiGrid.hpp"
#include "combigrid/combigrid/SerialCombiGrid.hpp"
#include "combigrid/fullgrid/CombiFullGrid.hpp"

#include "combigrid/combischeme/CombiSchemeBasis.hpp"
#include "combigrid/combischeme/CombiS_CT.hpp"

#include "combigrid/combischeme/CombiTS_CT.hpp"
#include "combigrid/combischeme/CombiArbitraryScheme.hpp"

#include "combigrid/domain/AbstractStretchingMaker.hpp"
#include "combigrid/domain/CombiGridDomain.hpp"
#include "combigrid/domain/CombiDomain1D.hpp"
#include "combigrid/domain/CombiTanStretching.hpp"
#include "combigrid/domain/CombiAtanSpecialStretching.hpp"
#include "combigrid/domain/CombiUniformStretching.hpp"

#include "combigrid/converter/CombiSGppConverter.hpp"

#endif /* COMBIGRID_HPP_ */
