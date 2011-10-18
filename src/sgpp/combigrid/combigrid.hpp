/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de

#ifndef COMBIGRID_HPP
#define COMBIGRID_HPP

#include "combigrid/combigrid/AbstractCombiGrid.hpp"
#include "combigrid/combigrid/SerialCombiGrid.hpp"
#include "combigrid/combigrid/AdaptiveSerialCombiGrid.hpp"
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



#endif /* COMBIGRID_HPP */
