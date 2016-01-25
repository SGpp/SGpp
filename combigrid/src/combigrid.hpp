// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_HPP_
#define COMBIGRID_HPP_

// --------- Core SGpp includes -----------
#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/base/grid/GridStorage.hpp"

// ----------- Combi Grid Includes ----------
#include "sgpp/combigrid/combigrid/AbstractCombiGrid.hpp"
#include "sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
#include "sgpp/combigrid/combigrid/AdaptiveSerialCombiGrid.hpp"
#include "sgpp/combigrid/combigrid/AdaptiveSerialCombiGridVariableCoefficients.hpp"
#include "sgpp/combigrid/fullgrid/CombiFullGrid.hpp"

#include "sgpp/combigrid/combischeme/CombiSchemeBasis.hpp"
#include "sgpp/combigrid/combischeme/CombiS_CT.hpp"

#include "sgpp/combigrid/combischeme/CombiTS_CT.hpp"
#include "sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp"

#include "sgpp/combigrid/domain/AbstractStretchingMaker.hpp"
#include "sgpp/combigrid/domain/CombiGridDomain.hpp"
#include "sgpp/combigrid/domain/CombiDomain1D.hpp"
#include "sgpp/combigrid/domain/CombiTanStretching.hpp"
#include "sgpp/combigrid/domain/CombiAtanSpecialStretching.hpp"
#include "sgpp/combigrid/domain/CombiUniformStretching.hpp"

#include "sgpp/combigrid/converter/CombiSGppConverter.hpp"

#include "sgpp/combigrid/operation/hash/OperationMatrixLTwoExplicitFullGrid.hpp"


#endif /* COMBIGRID_HPP_ */