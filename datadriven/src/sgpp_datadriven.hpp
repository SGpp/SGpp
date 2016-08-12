// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATADRIVEN_HPP
#define DATADRIVEN_HPP

#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostBase.hpp>
#include <sgpp/datadriven/algorithm/AlgorithmAdaBoostIdentity.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>
#include <sgpp/datadriven/algorithm/DMWeightMatrix.hpp>
#include <sgpp/datadriven/algorithm/DensitySystemMatrix.hpp>
#include <sgpp/datadriven/algorithm/test_dataset.hpp>

#include <sgpp/datadriven/functors/MultiGridRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/MultiSurplusRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/DataBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/GridPointBasedRefinementFunctor.hpp>
#include <sgpp/datadriven/functors/classification/ZeroCrossingRefinementFunctor.hpp>

#include <sgpp/datadriven/application/DensityEstimator.hpp>
#include <sgpp/datadriven/application/GaussianKDE.hpp>
#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/datadriven/application/LearnerSGDE.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRegularizationDiagonal.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTest.hpp>

#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimpleParameters.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformation.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationRosenblattTransformationKDE.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationTransformation1D.hpp>

#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMarginalize.hpp>

#include <sgpp/datadriven/tools/TypesDatadriven.hpp>

#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>

/* ************************
 * datamining
 * ************************/

// datasource
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceIterator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleDecorator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/SampleProvider.hpp>

#endif /* DATADRIVEN_HPP */
