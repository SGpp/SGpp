// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>
#include <sgpp/globaldef.hpp>

#if USE_OCL == 1
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
using sgpp::base::OCLOperationConfiguration;
#endif

#include <string>
#include <tuple>
#include <vector>

using sgpp::datadriven::BlacsProcessGrid;

std::string uncompressFile(std::string fileName);

sgpp::base::DataMatrix* readReferenceMatrix(sgpp::base::GridStorage& storage, std::string fileName);

void doRandomRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                         sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen,
                         sgpp::base::DataVector& alpha);

void doRandomRefinements(sgpp::base::AdaptivityConfiguration& adaptivityConfig,
                         sgpp::base::Grid& grid, sgpp::base::GridGenerator& gridGen);

double compareVectors(sgpp::base::DataVector& results, sgpp::base::DataVector& resultsCompare);

void compareDatasets(const std::vector<std::tuple<std::string, double>>& fileNamesError,
                     sgpp::base::GridType gridType, size_t level,
                     sgpp::datadriven::OperationMultipleEvalConfiguration configuration);

double compareToReference(sgpp::base::GridType gridType, const std::string& fileName, size_t level,
                          sgpp::datadriven::OperationMultipleEvalConfiguration configuration);

void compareDatasetsTranspose(const std::vector<std::tuple<std::string, double>>& fileNamesError,
                              sgpp::base::GridType gridType, size_t level,
                              sgpp::datadriven::OperationMultipleEvalConfiguration configuration);
double compareToReferenceTranspose(
    sgpp::base::GridType gridType, const std::string& fileName, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration);

void compareDatasetsDistributed(const std::vector<std::tuple<std::string, double>>& fileNamesError,
                                sgpp::base::GridType gridType, size_t level,
                                sgpp::datadriven::OperationMultipleEvalConfiguration configuration,
                                std::shared_ptr<BlacsProcessGrid> processGrid);

double compareToReferenceDistributed(
    sgpp::base::GridType gridType, const std::string& fileName, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration,
    std::shared_ptr<BlacsProcessGrid> processGrid);

void compareDatasetsTransposeDistributed(
    const std::vector<std::tuple<std::string, double>>& fileNamesError,
    sgpp::base::GridType gridType, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration,
    std::shared_ptr<BlacsProcessGrid> processGrid);

double compareToReferenceTransposeDistributed(
    sgpp::base::GridType gridType, const std::string& fileName, size_t level,
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration,
    std::shared_ptr<BlacsProcessGrid> processGrid);

#if USE_OCL == 1

std::shared_ptr<OCLOperationConfiguration> getConfigurationDefaultsSingleDevice();

std::shared_ptr<OCLOperationConfiguration> getConfigurationDefaultsMultiDevice();

std::shared_ptr<OCLOperationConfiguration> getConfigurationDefaultsMultiPlatform();
#endif
