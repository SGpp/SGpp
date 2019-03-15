/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * Created by Bountos Nikolaos on 12/14/18
 */

#include <sgpp/datadriven/datamining/modules/fitting/PDFCombigrid.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>


PDFCombigrid::PDFCombigrid(const sgpp::datadriven::FitterConfigurationDensityEstimation &conf) {
    ModelFittingDensityEstimationOnOff();
    this->config = std::unique_ptr<sgpp::datadriven::FitterConfiguration>(
            std::make_unique<sgpp::datadriven::FitterConfigurationDensityEstimation>(conf));
    //dimensions = dataset->getDimension()+1;
    //model = new sgpp::datadriven::PDFFitter();
    //model->setDataset(dataset);
    //model->setConfiguration(*miner->getModel());
    this->level = conf.getGridConfig().level_ - 1 ;
    this->numthreads = conf.getGridConfig().threads;
    if ( this->numthreads > 1 )
        this->parallel = true;
    else
        this->parallel = false;
}

PDFCombigrid::PDFCombigrid(int lev, int threads, std::string conf, bool par) {
    DensityEstimationMinerFactory factory;
    configuration = conf;
    miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner(conf));
    sgpp::datadriven::DataSourceSplitting* datasource;
    sgpp::datadriven::DataSourceCrossValidation* dcross;
    sgpp::datadriven::DataMiningConfigParser  parser(conf);
    sgpp::datadriven::DataSourceConfig dconfig{};
    parser.getDataSourceConfig(dconfig, dconfig);
    dataset = new sgpp::datadriven::Dataset();
    sgpp::datadriven::DataSourceBuilder builder;
    sgpp::datadriven::DataSourceBuilder b = builder.withPath(dconfig.filePath);
    if (parser.hasFitterConfigCrossValidation()) {
        CrossvalidationConfiguration crossValidationconfig{};
        parser.getFitterCrossvalidationConfig(crossValidationconfig, crossValidationconfig);
        dcross = b.crossValidationFromConfig(dconfig, crossValidationconfig);
        dataset = dcross->getAllSamples();
    } else {
        datasource = b.splittingFromConfig(dconfig);
        dataset = datasource->getAllSamples();
    }
    dimensions = dataset->getDimension()+1;
    model = new sgpp::datadriven::PDFFitter();
    model->setDataset(dataset);
    model->setConfiguration(*miner->getModel());
    this->parallel = par;
    this->level = lev;
    this->numthreads = threads;
}

void PDFCombigrid::fit() {
    // grid function
    dimensions = dataset->getDimension()+1;
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
        /**
        * The miner object is constructed by the factory from a supplied configuration file.
        */
        auto parallel = this->parallel;
        sgpp::datadriven::PDFFitter *modelp = nullptr;
        if (parallel) {
            modelp = new sgpp::datadriven::PDFFitter();
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
            model->fit(*dataset, grids, an, true);
        } else {
            modelp->fit(*dataset, grids, an, true);
        }

        while (it.isValid()) {
            auto ik = it.getMultiIndex();
            // evaluate grid point
            double value = 0.0;
            if (!parallel) {
                value = model->evaluate(grid->getGridPoint(ik));
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

    /**
     * We have to specify if the function always produces the same value for the same grid points.
     * This can make the storage smaller if the grid points are nested. In this implementation, this
     * is true. However, it would be false in the PDE case, so we set it to false here.
     */
    bool exploitNesting = false;


    operation = std::make_shared<sgpp::combigrid::CombigridOperation>(
            grids, evaluators, levelManager, gf, exploitNesting);
}

void PDFCombigrid::update(sgpp::datadriven::Dataset &newDataset) {
    std::cout<<"PAME GERA"<<std::endl;
    dataset = &newDataset;
    std::cout<<" Dimension "<< dataset->getDimension()<<std::endl;
    fit();
    std::cout<< " Out of fit "<<std::endl;
}

double PDFCombigrid::evaluate(std::vector<double> test_points) {
    sgpp::base::DataVector parameter(test_points);
    double result;
    if (!parallel)
        result = operation->evaluate(this->level, parameter);
    else
        result = operation->evaluateParallel(this->level, parameter, numthreads);

    return result;
}

void PDFCombigrid::evaluate(DataMatrix &samples, DataVector &results) {
    for ( auto i=0; i <= samples.getNrows(); i++)
    {
        std::vector<double> row;
        size_t current = 0;
        for ( auto j= i * samples.getNcols(); j <= (i + 1)*samples.getNcols(); j++)
        {
            row.push_back(samples[j]);
        }
        results.push_back(evaluate(row));
    }
}