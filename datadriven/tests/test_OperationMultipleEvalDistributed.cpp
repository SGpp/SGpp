// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_SCALAPACK

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalScalapack/OperationMultipleEvalDistributed.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <string>
#include <tuple>
#include <vector>

#include "test_datadrivenCommon.hpp"

using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::OperationMultipleEvalConfiguration;

namespace TestOperationMultipleEvalDistributedFixture {
struct FilesNamesAndErrorFixture {
  FilesNamesAndErrorFixture() {}
  ~FilesNamesAndErrorFixture() {}

  std::vector<std::tuple<std::string, double>> fileNamesErrorDouble = {
      std::tuple<std::string, double>("datadriven/datasets/friedman/friedman2_4d_10000.arff.gz",
                                      1E-15),
      std::tuple<std::string, double>("datadriven/datasets/friedman/friedman1_10d_2000.arff.gz",
                                      1E-15)};

  uint32_t level = 5;
};
}  // namespace TestOperationMultipleEvalDistributedFixture

BOOST_FIXTURE_TEST_SUITE(testOperationMultipleEvalDistributed,
                         TestOperationMultipleEvalDistributedFixture::FilesNamesAndErrorFixture)
#ifdef ZLIB

BOOST_AUTO_TEST_CASE(multLinear) {
  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::SCALAPACK);

  std::shared_ptr<BlacsProcessGrid> processGrid = std::make_shared<BlacsProcessGrid>();

  compareDatasetsDistributed(fileNamesErrorDouble, sgpp::base::GridType::Linear, level,
                             configuration, processGrid);
}

BOOST_AUTO_TEST_CASE(multModLinear) {
  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::SCALAPACK);

  std::shared_ptr<BlacsProcessGrid> processGrid = std::make_shared<BlacsProcessGrid>();

  compareDatasetsDistributed(fileNamesErrorDouble, sgpp::base::GridType::ModLinear, level,
                             configuration, processGrid);
}

BOOST_AUTO_TEST_CASE(multTransposeLinear) {
  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::SCALAPACK);

  std::shared_ptr<BlacsProcessGrid> processGrid = std::make_shared<BlacsProcessGrid>();

  compareDatasetsTransposeDistributed(fileNamesErrorDouble, sgpp::base::GridType::Linear, level,
                                      configuration, processGrid);
}

BOOST_AUTO_TEST_CASE(multTransposeModLinear) {
  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::SCALAPACK);

  std::shared_ptr<BlacsProcessGrid> processGrid = std::make_shared<BlacsProcessGrid>();

  compareDatasetsTransposeDistributed(fileNamesErrorDouble, sgpp::base::GridType::ModLinear, level,
                                      configuration, processGrid);
}

#else

BOOST_AUTO_TEST_CASE(noZlibOpMultipleEvalDistributed) {
  BOOST_TEST_MESSAGE(
      "Zlib not found, cannot load files for test of OperationMultipleEvalDistributed");
}
#endif /* ZLIB */

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_SCALAPACK */
