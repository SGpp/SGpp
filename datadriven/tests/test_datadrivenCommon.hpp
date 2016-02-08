/*
 * test_datadrivenCommon.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#if USE_OCL==1
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#endif
#include <sgpp/globaldef.hpp>

std::string uncompressFile(std::string fileName);

SGPP::base::DataMatrix* readReferenceMatrix(SGPP::base::GridStorage* storage,
    std::string fileName);

void doRandomRefinements(SGPP::base::AdpativityConfiguration& adaptConfig,
                         SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen,
                         SGPP::base::DataVector& alpha);

void doRandomRefinements(SGPP::base::AdpativityConfiguration& adaptConfig,
                         SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen);

double compareVectors(SGPP::base::DataVector& results,
                      SGPP::base::DataVector& resultsCompare);

double compareToReference(SGPP::base::GridType gridType, std::string fileName,
                          size_t level,
                          SGPP::datadriven::OperationMultipleEvalConfiguration configuration,
                          size_t numRefinements = 1);

double compareToReferenceTranspose(SGPP::base::GridType gridType,
                                   std::string fileName, size_t level,
                                   SGPP::datadriven::OperationMultipleEvalConfiguration configuration);
#if USE_OCL==1
SGPP::base::OCLOperationConfiguration getConfigurationDefaultsSingleDevice();

SGPP::base::OCLOperationConfiguration getConfigurationDefaultsMultiDevice();

SGPP::base::OCLOperationConfiguration getConfigurationDefaultsMultiPlatform();
#endif
