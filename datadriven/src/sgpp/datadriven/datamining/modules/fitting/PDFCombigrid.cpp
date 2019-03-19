/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Created by Bountos Nikolaos on 12/14/18
 */

#include <sgpp/datadriven/datamining/modules/fitting/PDFCombigrid.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <string>
#include <vector>

std::mutex mymut;

PDFCombigrid::PDFCombigrid(const sgpp::datadriven::FitterConfigurationDensityEstimation &conf) {
    ModelFittingDensityEstimationOnOff();
    this->config = std::unique_ptr<sgpp::datadriven::FitterConfiguration>(
            std::make_unique<sgpp::datadriven::FitterConfigurationDensityEstimation>(conf));
    this->level = conf.getGridConfig().level_ - 1;
    this->numthreads = conf.getGridConfig().threads;
    if ( this->numthreads > 1 )
        this->parallel = true;
    else
        this->parallel = false;
    current = 0;
    fitted = false;
}

void PDFCombigrid::fit() {
    // grid function
    dimensions = dataset->getDimension() + 1;
    model = new sgpp::datadriven::PDFFitter();
    model->setDataset(dataset);
    model->setConfiguration(*this);
    sgpp::combigrid::GridFunction gf([this](std::shared_ptr<sgpp::combigrid::TensorGrid> grid) {
        // We store the results for each grid point, encoded by a MultiIndex, in a TreeStorage
        auto result = std::make_shared<sgpp::combigrid::TreeStorage<double> >(dimensions);
        // Creates an iterator that yields all multi-indices of grid points in the grid.
        sgpp::combigrid::MultiIndexIterator it(grid->numPoints());
        auto  num = grid->numPoints();
        // create grid based on datavector a
        // make a model out of it
        // fit it
        // evaluate
        auto parallel = this->parallel;
        sgpp::datadriven::PDFFitter *modelp = nullptr;
        if (!fitted) {
            mymut.lock();
            modelp = modelpool[current];
            current++;
            mymut.unlock();
            modelp->setDataset(this->dataset);
            modelp->setConfiguration(*model);
        }
        std::unique_ptr<sgpp::base::Grid> grids(sgpp::base::Grid::createLinearGrid(num.size()));
        auto an = grid->getLevel();
        std::string m;
        // increase level for datadriven module consistency
        for (size_t  i = 0; i < an.size(); i++) {
            an[i]++;
            m = m+""+std::to_string(an[i]);
        }
        // fit dataset to Grid
        if (!parallel) {
            if (!fitted) {
                modelp->fit(*dataset, grids, an, true);
                models.insert(std::make_pair(m, modelp));
            } else {
                modelp = models.at(m);
            }
        } else {
            if (!fitted) {
                modelp->fit(*dataset, grids, an, true);
                mymut.lock();
                models.insert(std::make_pair(m, modelp));
                mymut.unlock();
            } else {
                modelp = models.at(m);
            }
        }
        while (it.isValid()) {
            auto ik = it.getMultiIndex();
            // evaluate grid point
            double value = 0.0;
            if (!parallel) {
                value = modelp->evaluate(grid->getGridPoint(ik));
            } else {
                value = modelp->evaluate(grid->getGridPoint(ik));
            }
            // Store the result at the multi index encoding the grid point
            result->set(it.getMultiIndex(), value);
            it.moveToNext();
        }
        return result;
    });

    sgpp::combigrid::CombiHierarchies::Collection grids(
            dimensions, sgpp::combigrid::CombiHierarchies::expUniform());
    sgpp::combigrid::CombiEvaluators::Collection evaluators(
            dimensions, sgpp::combigrid::CombiEvaluators::linearInterpolation());
    std::shared_ptr<sgpp::combigrid::LevelManager> levelManager(
            new sgpp::combigrid::AveragingLevelManager());

    bool exploitNesting = false;

    operation = std::make_shared<sgpp::combigrid::CombigridOperation>(
            grids, evaluators, levelManager, gf, exploitNesting);
}

void PDFCombigrid::update(sgpp::datadriven::Dataset &newDataset) {
    dataset = &newDataset;
    fit();
    numcomponents = pow(2, pow(this->level, this->dimensions));
    for (int i = 0 ; i < numcomponents; i++) {
        sgpp::datadriven::PDFFitter *modelp =  new sgpp::datadriven::PDFFitter();
        modelpool.push_back(modelp);
    }
    auto& gridConfig = this->config->getGridConfig();
    grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};
}

double PDFCombigrid::evaluate(std::vector<double> test_points) {
    sgpp::base::DataVector parameter(test_points);
    double result;
    if (!parallel)
        result = operation->evaluate(this->level, parameter);
    else
        result = operation->evaluateParallel(this->level, parameter, numthreads);
    fitted = true;
    return result;
}

void PDFCombigrid::evaluate(DataMatrix &samples, DataVector &results) {
    for ( auto i = 0; i < samples.getNrows(); i++) {
        std::vector<double> row;
        for (auto j=0; j <= samples.getNcols() ; j++) {
            auto index = i*(samples.getNcols()+1)+j;
            row.push_back(samples[index]);
        }
        DataVector param(row);
        results[i]=(evaluate(param));
    }
}

bool PDFCombigrid::refine() {
    return false;
}
