/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Created by Bountos Nikolaos on 12/14/18
 */

#pragma once

#ifndef SGPP_PDECOMBIGRID_H
#define SGPP_PDECOMBIGRID_H



#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/DensityEstimationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/grid/TensorGrid.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/CombigridEvaluator.hpp>
#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/LinearInterpolationEvaluator.hpp>
#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/GeneralFunction.cpp>
#include <sgpp/datadriven/datamining/modules/fitting/PDFFitter.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>


using sgpp::datadriven::DensityEstimationMinerFactory;
using sgpp::datadriven::SparseGridMiner;
using sgpp::combigrid::MultiIndex;
using sgpp::datadriven::ModelFittingDensityEstimationOnOff;

class PDFCombigrid : public ModelFittingDensityEstimationOnOff {
 private:
    size_t dimensions;
    sgpp::datadriven::Dataset* dataset;
    std::unique_ptr<SparseGridMiner> miner;
    sgpp::datadriven::PDFFitter* model;
    std::string configuration;
    bool parallel;
    int level;
    int numthreads;

    std::shared_ptr<sgpp::combigrid::CombigridOperation> operation;
    bool fitted;

public:
    PDFCombigrid(int lev, int threads, std::string conf, bool par);
    PDFCombigrid(const sgpp::datadriven::FitterConfigurationDensityEstimation& conf);
    void update(sgpp::datadriven::Dataset& newDataset);
    void fit();
    double evaluate(std::vector<double> test_points);
    /**
     * Evaluate the fitted density on a set of data points - requires a trained grid.
     * @param samples matrix where each row represents a sample and the columns contain the
     * coordinates in all dimensions of that sample.
     * @param results vector where each row will contain the evaluation of the respective sample on
     * the current model.
        */
    void evaluate(DataMatrix& samples, DataVector& results) override;
};


#endif  // SGPP_PDECOMBIGRID_H
