/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Created by Bountos Nikolaos on 12/14/18
 */
#include <sgpp/datadriven/datamining/modules/fitting/PDFCombigrid.hpp>
#include<vector>
#include <string>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/UniversalMinerFactory.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
int main(int argc, char **argv) {
    std::string s = "../../../../myconfig.json";

    sgpp::datadriven::UniversalMinerFactory factory;

    auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner("../../../../myconfig.json"));
    std::cout<< typeid(miner->getModel()).name()<<std::endl;
    miner->learn(false);
    std::vector<double>  a;
    a.push_back(0.1);
    a.push_back(0.2);
    DataVector param(a);
    std::cout<<std::endl<<std::endl<<miner->getModel()->evaluate(param);
   // auto miner2 = std::unique_ptr<SparseGridMiner>(factory.buildMiner("../../../../miner.json"));
    //miner2->learn(false);
    //std::cout<<std::endl<<std::endl<<miner2->getModel()->evaluate(param);

}
