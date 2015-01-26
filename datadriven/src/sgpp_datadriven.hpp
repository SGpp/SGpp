/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de

#ifndef DATADRIVEN_HPP
#define DATADRIVEN_HPP


#include "sgpp/datadriven/algorithm/test_dataset.hpp"
#include "sgpp/datadriven/algorithm/DMSystemMatrix.hpp"
#include "sgpp/datadriven/algorithm/DMWeightMatrix.hpp"
#include "sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp"
#include "sgpp/datadriven/algorithm/AlgorithmAdaBoostIdentity.hpp"
#include "sgpp/datadriven/algorithm/DensitySystemMatrix.hpp"

#include "sgpp/datadriven/application/Learner.hpp"
#include "sgpp/datadriven/application/LearnerDensityBased.hpp"
#include "sgpp/datadriven/application/LearnerDensityBasedReg.hpp"
#include "sgpp/datadriven/application/LearnerSGD.hpp"
#include "sgpp/datadriven/application/LearnerOnlineSGD.hpp"
#include "sgpp/datadriven/application/LearnerDensityCluster.hpp"

#include "sgpp/datadriven/operation/OperationRegularizationDiagonal.hpp"
#include "sgpp/datadriven/operation/OperationTest.hpp"

#include "sgpp/datadriven/tools/ARFFTools.hpp"

#include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp"
#include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/CommonParameters.hpp"
#include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp"
#include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimpleParameters.hpp"

#include "sgpp/datadriven/tools/TypesDatadriven.hpp"


#endif /* DATADRIVEN_HPP */
