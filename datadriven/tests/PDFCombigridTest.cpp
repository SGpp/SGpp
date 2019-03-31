/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Created by Bountos Nikolaos on 12/14/18
 */
#if defined USE_SGDECOMBI && defined USE_GSL
#include <sgpp/datadriven/datamining/modules/fitting/PDFCombigrid.hpp>
#include <sgpp/datadriven/datamining/builder/UniversalMinerFactory.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/DensityEstimationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>

using sgpp::datadriven::DensityEstimationMinerFactory;
using sgpp::datadriven::SparseGridMiner;
using sgpp::datadriven::ModelFittingBase;
using sgpp::datadriven::CSVFileSampleProvider;
using sgpp::datadriven::DataVector;

double testParallel(int dim) {
    sgpp::datadriven::UniversalMinerFactory factory;
    auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner("PDFCombi.json"));
    miner->learn(false);
    auto minertrue = std::unique_ptr<SparseGridMiner>(factory.buildMiner("miner.json"));
    minertrue->learn(false);
    std::vector<double> a;
    for (auto i = 0; i < dim ; i++) {
        a.push_back(0.1);
    }
    std::cout << " Result  Predicted " << dynamic_cast<PDFCombigrid*>
    ((miner->getModel()))->evaluate(DataVector(a))
    << " Ground truth " << minertrue->getModel()->evaluate(DataVector(a));
    return dynamic_cast<PDFCombigrid*>((miner->getModel()))->
    evaluate(DataVector(a)) - minertrue->getModel()->evaluate(DataVector(a));
}

double testSequential(int dim) {
    sgpp::datadriven::UniversalMinerFactory factory;
    auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner("PDFSequential.json"));
    miner->learn(false);
    auto minertrue = std::unique_ptr<SparseGridMiner>(factory.buildMiner("miner.json"));
    minertrue->learn(false);
    std::vector<double> a;
    for (auto i = 0; i < dim ; i++) {
        a.push_back(0.1);
    }
    std::cout << " Result  Predicted " << dynamic_cast<PDFCombigrid*>
    ((miner->getModel()))->evaluate(DataVector(a))
    << " Ground truth " << minertrue->getModel()->evaluate(DataVector(a));
    return dynamic_cast<PDFCombigrid*>((miner->getModel()))->
    evaluate(DataVector(a)) - minertrue->getModel()->evaluate(DataVector(a));
}

BOOST_AUTO_TEST_SUITE(testPDFCombigridTest)

BOOST_AUTO_TEST_CASE(testResult) {
  double diff = testParallel(2);
  BOOST_CHECK(diff < 2);
}

BOOST_AUTO_TEST_CASE(testResultSequential) {
  double diff = testSequential(2);
  BOOST_CHECK(diff < 2);
}


BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_SGDECOMBI */

