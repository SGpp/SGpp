// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef DISABLED

#if USE_OCL == 1

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <zlib.h>

#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "test_datadrivenCommon.hpp"
#include "sgpp/globaldef.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/base/tools/ConfigurationParameters.hpp"

BOOST_AUTO_TEST_SUITE(TestStreamingOCLMultTranspose)

BOOST_AUTO_TEST_CASE(Simple) {
  std::vector<std::tuple<std::string, double> > fileNamesError = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17)};

  uint32_t level = 5;

  SGPP::base::OCLOperationConfiguration parameters = getConfigurationDefaultsSingleDevice();
  parameters.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
  parameters.replaceIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
  parameters.replaceIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 1ul);
  parameters.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", 1ul);
  parameters.replaceTextAttr("KERNEL_STORE_DATA", "array");

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceTranspose(SGPP::base::GridType::Linear,
                                             std::get<0>(fileNameError), level, configuration);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
  }
}

BOOST_AUTO_TEST_CASE(Blocking) {
  std::vector<std::tuple<std::string, double> > fileNamesError = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17)};

  uint32_t level = 6;

  SGPP::base::OCLOperationConfiguration parameters = getConfigurationDefaultsSingleDevice();
  parameters.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
  parameters.replaceIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
  parameters.replaceIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 4ul);
  parameters.replaceTextAttr("KERNEL_STORE_DATA", "register");
  parameters.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceTranspose(SGPP::base::GridType::Linear,
                                             std::get<0>(fileNameError), level, configuration);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
  }
}

BOOST_AUTO_TEST_CASE(MultiDevice) {
  //  std::string fileName = "friedman2_90000.arff";
  //    std::string fileName = "debugging.arff";
  std::vector<std::tuple<std::string, double> > fileNamesError = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman_4d.arff.gz", 1E-13),
      std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 1E-17)};

  uint32_t level = 6;

  SGPP::base::OCLOperationConfiguration parameters = getConfigurationDefaultsMultiDevice();
  parameters.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);
  parameters.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
  parameters.replaceIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
  parameters.replaceIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 4ul);
  parameters.replaceTextAttr("KERNEL_STORE_DATA", "register");

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceTranspose(SGPP::base::GridType::Linear,
                                             std::get<0>(fileNameError), level, configuration);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
  }
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {
  std::vector<std::tuple<std::string, double> > fileNamesError = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman_4d.arff.gz", 1E+4),
      std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 20.0)};

  uint32_t level = 6;

  SGPP::base::OCLOperationConfiguration parameters = getConfigurationDefaultsSingleDevice();
  parameters.replaceTextAttr("INTERNAL_PRECISION", "float");
  parameters.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
  parameters.replaceIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
  parameters.replaceIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 1ul);
  parameters.replaceTextAttr("KERNEL_STORE_DATA", "array");
  parameters.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceTranspose(SGPP::base::GridType::Linear,
                                             std::get<0>(fileNameError), level, configuration);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
  }
}

BOOST_AUTO_TEST_CASE(BlockingSinglePrecision) {
  std::vector<std::tuple<std::string, double> > fileNamesError = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman_4d.arff.gz", 1E+4),
      std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 20.0)};

  uint32_t level = 6;

  SGPP::base::OCLOperationConfiguration parameters = getConfigurationDefaultsSingleDevice();
  parameters.replaceTextAttr("INTERNAL_PRECISION", "float");
  parameters.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
  parameters.replaceIDAttr("KERNEL_DATA_BLOCKING_SIZE", 4ul);
  parameters.replaceIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 1ul);
  parameters.replaceTextAttr("KERNEL_STORE_DATA", "array");
  parameters.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceTranspose(SGPP::base::GridType::Linear,
                                             std::get<0>(fileNameError), level, configuration);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
  }
}

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {
  std::vector<std::tuple<std::string, double> > fileNamesError = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman_4d.arff.gz", 1E+4),
      std::tuple<std::string, double>("datadriven/tests/data/friedman_10d.arff.gz", 20.0)};

  uint32_t level = 6;

  SGPP::base::OCLOperationConfiguration parameters = getConfigurationDefaultsMultiDevice();
  parameters.replaceTextAttr("INTERNAL_PRECISION", "float");
  parameters.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
  parameters.replaceIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
  parameters.replaceIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 1ul);
  parameters.replaceTextAttr("KERNEL_STORE_DATA", "array");
  parameters.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCL, parameters);

  for (std::tuple<std::string, double> fileNameError : fileNamesError) {
    double mse = compareToReferenceTranspose(SGPP::base::GridType::Linear,
                                             std::get<0>(fileNameError), level, configuration);
    BOOST_CHECK(mse < std::get<1>(fileNameError));
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
