// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB
#if USE_OCL == 1

#define BOOST_TEST_DYN_LINK
#include <zlib.h>
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/base/tools/ConfigurationParameters.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/globaldef.hpp"
#include "test_datadrivenCommon.hpp"

using sgpp::datadriven::OperationMultipleEvalConfiguration;
using sgpp::datadriven::OperationMultipleEvalType;
using sgpp::datadriven::OperationMultipleEvalSubType;

namespace TestStreamingModOCLMaskMultiPlatformMultTransposeFixture {
struct FilesNamesAndErrorFixture {
  FilesNamesAndErrorFixture() {}
  ~FilesNamesAndErrorFixture() {}

  std::vector<std::tuple<std::string, double>> fileNamesErrorDouble = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman2_4d_10000.arff.gz", 1E-16),
      std::tuple<std::string, double>("datadriven/tests/data/friedman1_10d_2000.arff.gz", 1E-20)};

  std::vector<std::tuple<std::string, double>> fileNamesErrorFloat = {
      std::tuple<std::string, double>("datadriven/tests/data/friedman2_4d_10000.arff.gz", 1E+5),
      std::tuple<std::string, double>("datadriven/tests/data/friedman1_10d_2000.arff.gz", 1E-3)};

  uint32_t level = 5;
};
}  // namespace TestStreamingModOCLMaskMultiPlatformMultTransposeFixture

BOOST_FIXTURE_TEST_SUITE(
    TestStreamingModOCLMaskMultiPlatformMultTranspose,
    TestStreamingModOCLMaskMultiPlatformMultTransposeFixture::FilesNamesAndErrorFixture)

BOOST_AUTO_TEST_CASE(Simple) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();

  (*parameters)["INTERNAL_PRECISION"].set("double");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(1));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorDouble, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(Local) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();

  (*parameters)["INTERNAL_PRECISION"].set("double");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorDouble, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(Blocking) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();

  (*parameters)["INTERNAL_PRECISION"].set("double");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
    deviceNode.replaceIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(2));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorDouble, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(MultiDevice) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsMultiDevice();

  (*parameters)["INTERNAL_PRECISION"].set("double");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
    deviceNode.replaceIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(2));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorDouble, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(MultiPlatform) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsMultiPlatform();

  (*parameters)["INTERNAL_PRECISION"].set("double");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
    deviceNode.replaceIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(2));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorDouble, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(SimpleSinglePrecision) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();

  (*parameters)["INTERNAL_PRECISION"].set("float");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(1));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorFloat, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(LocalSinglePrecision) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();

  (*parameters)["INTERNAL_PRECISION"].set("float");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorFloat, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(BlockingSinglePrecision) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsSingleDevice();

  (*parameters)["INTERNAL_PRECISION"].set("float");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
    deviceNode.replaceIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(2));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorFloat, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(MultiDeviceSinglePrecision) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsMultiDevice();

  (*parameters)["INTERNAL_PRECISION"].set("float");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
    deviceNode.replaceIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(2));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorFloat, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_CASE(MultiPlatformSinglePrecision) {
  std::shared_ptr<sgpp::base::OCLOperationConfiguration> parameters =
      getConfigurationDefaultsMultiPlatform();

  (*parameters)["INTERNAL_PRECISION"].set("float");

  std::vector<std::reference_wrapper<json::Node>> deviceNodes = parameters->getAllDeviceNodes();
  for (json::Node &deviceNode : deviceNodes) {
    deviceNode.replaceIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
    deviceNode.replaceIDAttr("KERNEL_MAX_DIM_UNROLL", UINT64_C(10));
    deviceNode.replaceIDAttr("KERNEL_TRANS_GRID_BLOCK_SIZE", UINT64_C(2));
    deviceNode.replaceIDAttr("KERNEL_STORE_DATA", "register");
  }

  OperationMultipleEvalConfiguration configuration(
      OperationMultipleEvalType::STREAMING, OperationMultipleEvalSubType::OCLMASKMP, *parameters);

  compareDatasetsTranspose(fileNamesErrorFloat, sgpp::base::GridType::ModLinear, level,
                           configuration);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
