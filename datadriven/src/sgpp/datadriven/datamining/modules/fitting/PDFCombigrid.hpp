// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once


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
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDE.hpp>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>


using sgpp::datadriven::DensityEstimationMinerFactory;
using sgpp::datadriven::SparseGridMiner;
using sgpp::combigrid::MultiIndex;
using sgpp::datadriven::ModelFittingDensityEstimationOnOff;

class PDFCombigrid : public ModelFittingDensityEstimationOnOff {
 private:
    // number of dimensions
    size_t dimensions;
    // Model's dataset
    sgpp::datadriven::Dataset* dataset;
    // PDFFitter main model
    sgpp::datadriven::PDFFitter* model;
    // model's configuration
    std::string configuration;
    // Parallel option
    bool parallel;
    // Grid's Level
    int level;
    // number of threads
    int numthreads;
    // Combigrid's Operation object
    #ifdef USE_SGDECOMBI
    std::shared_ptr<sgpp::combigrid::CombigridOperation> operation;
    #endif
    // True if model is in evaluation mode, False if the model has not been fitted
    bool fitted;
    // Hash map containing the model for each subspace
    std::unordered_map<std::string, sgpp::datadriven::PDFFitter*> models;
    // Initialized  model pool to be used by each subspace.
    std::vector<sgpp::datadriven::PDFFitter*> modelpool;
    // Estimated number of components
    int numcomponents;
    // Index of component subspace used to grab model from pool.
    int current;

 public:
    explicit PDFCombigrid(const sgpp::datadriven::FitterConfigurationDensityEstimation& conf);
    /**
     * Destructor to clear the pool of models
     */
    ~PDFCombigrid();
    /**
     * Model update on the given dataset
     * @param newDataset
     */
    void update(sgpp::datadriven::Dataset& newDataset);
    /**
     * Fit themodel
     */
    void fit();
    /**
     * Evaluate function on test_points
     * @param test_points
     * @return evaluation
     */
    double evaluate(std::vector<double> test_points);
    /**
     * Function to skip refine in scorer
     * @return false
     */
    bool refine() override;
    /**
     * Evaluate the fitted density on a set of data points - requires a trained grid.
     * @param samples matrix where each row represents a sample and the columns contain the
     * coordinates in all dimensions of that sample.
     * @param results vector where each row will contain the evaluation of the respective sample on
     * the current model.
        */
    void evaluate(DataMatrix& samples, DataVector& results) override;
};
