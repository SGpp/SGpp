// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATADRIVEN_HPP
#define DATADRIVEN_HPP

#include <sgpp/datadriven/algorithm/test_dataset.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>
#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>
#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp>
#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostIdentity.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/algorithm/DBMatOffline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnline.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <sgpp/datadriven/algorithm/DBMatDecompMatrixSolver.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSBackSub.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSEigen.hpp>
#include <sgpp/datadriven/algorithm/DBMatDMSChol.hpp>
#include <sgpp/datadriven/algorithm/DBMatDensityConfiguration.hpp> 
#include <sgpp/datadriven/algorithm/ConvergenceMonitor.hpp>

#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>


#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/datadriven/application/GaussianKDE.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>
#include <sgpp/datadriven/application/LearnerSGDEOnOff.hpp>
#include <sgpp/datadriven/application/PrimalDualSVM.hpp>
#include <sgpp/datadriven/application/LearnerSVM.hpp>
#include <sgpp/datadriven/application/LearnerSGD.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRegularizationDiagonal.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>

#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimpleParameters.hpp>
#include <sgpp/datadriven/operation/hash/OperationEvalSGKernel/OperationEvalSGKernel.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationTransformation1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>

#include <sgpp/datadriven/tools/TypesDatadriven.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

#endif /* DATADRIVEN_HPP */
