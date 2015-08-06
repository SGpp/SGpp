// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
#include "sgpp/datadriven/application/DensityEstimator.hpp"
#include "sgpp/datadriven/application/GaussianKDE.hpp"
#include "sgpp/datadriven/application/LearnerSGDE.hpp"

#include "sgpp/datadriven/operation/hash/OperationRegularizationDiagonal.hpp"
#include "sgpp/datadriven/operation/hash/OperationTest.hpp"

#include "sgpp/datadriven/tools/ARFFTools.hpp"

#ifdef __AVX__
#include "sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp"
#include "sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp"
#include "sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp"
#include "sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimpleParameters.hpp"
#endif

#include "sgpp/datadriven/tools/TypesDatadriven.hpp"

#include "sgpp/datadriven/DatadrivenOpFactory.hpp"


#endif /* DATADRIVEN_HPP */
