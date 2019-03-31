// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#include <sgpp/datadriven/datamining/modules/fitting/PDFCombigrid.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/UniversalMinerFactory.hpp>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
int main(int argc, char **argv) {
    sgpp::datadriven::UniversalMinerFactory factory;
    auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner("./PDFCombi.json"));
    miner->learn(false);
}
