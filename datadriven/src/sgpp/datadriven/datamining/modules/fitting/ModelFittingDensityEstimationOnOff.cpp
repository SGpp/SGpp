/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ModelFittingDensityEstimation.cpp
 *
 * Created on: Jan 02, 2018
 *     Author: Kilian Röhner
 */

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatDatabase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOff.hpp>
#include <string>
#include <vector>
#include <list>

#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingClassification.hpp>


using sgpp::base::Grid;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::SurplusRefinementFunctor;
using sgpp::base::RefinementFunctor;
using sgpp::base::SurplusVolumeRefinementFunctor;
using sgpp::base::RefinementFunctorType;

using sgpp::base::application_exception;

namespace sgpp {
namespace datadriven {

ModelFittingDensityEstimationOnOff::ModelFittingDensityEstimationOnOff(
    const FitterConfigurationDensityEstimation& config) : ModelFittingDensityEstimation() {
  this->config = std::unique_ptr<FitterConfiguration>(
      std::make_unique<FitterConfigurationDensityEstimation>(config));
}

// TODO(lettrich): exceptions have to be thrown if not valid.
double ModelFittingDensityEstimationOnOff::evaluate(const DataVector& sample) {
  return online->eval(alpha, sample, *grid);
}

// TODO(lettrich): exceptions have to be thrown if not valid.
void ModelFittingDensityEstimationOnOff::evaluate(DataMatrix& samples, DataVector& results) {
    online->eval(alpha, samples, results, *grid);
}

void ModelFittingDensityEstimationOnOff::fit(Dataset& newDataset) {
  dataset = &newDataset;
  fit(newDataset.getData());
}

void ModelFittingDensityEstimationOnOff::fit(DataMatrix& newDataset) {
  // Get configurations
  auto& databaseConfig = this->config->getDatabaseConfig();
  auto& gridConfig = this->config->getGridConfig();
  auto& refinementConfig = this->config->getRefinementConfig();
  auto& regularizationConfig = this->config->getRegularizationConfig();
  auto& densityEstimationConfig = this->config->getDensityEstimationConfig();
  auto& geometryConfig = this->config->getGeometryConfig();

  // clear model
  reset();

  // build grid
  gridConfig.dim_ = newDataset.getNcols();
  std::cout << "Dataset dimension " << gridConfig.dim_ << std::endl;
  // TODO(fuchsgruber): Support for geometry aware sparse grids (pass interactions from config?)
  // grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};
  // grid = std::unique_ptr<Grid>{buildGrid(gridConfig)};
  grid = std::unique_ptr<Grid>{buildGrid(gridConfig, geometryConfig)};

  // build surplus vector
  alpha = DataVector{grid->getSize()};

  std::cout << grid->getSize() << std::endl;

  // Build the offline instance first
  DBMatOffline *offline = nullptr;

  // Intialize database if it is provided
  if (!databaseConfig.filepath.empty()) {
    datadriven::DBMatDatabase database(databaseConfig.filepath);
    // Check if database holds a fitting lhs matrix decomposition
    if (database.hasDataMatrix(gridConfig, refinementConfig, regularizationConfig,
        densityEstimationConfig)) {
      std::string offlineFilepath = database.getDataMatrix(gridConfig, refinementConfig,
          regularizationConfig, densityEstimationConfig);
      offline = DBMatOfflineFactory::buildFromFile(offlineFilepath);
    }
  }

  // Build and decompose offline object if not loaded from database
  if (offline == nullptr) {
    // Build offline object by factory, build matrix and decompose
    offline = DBMatOfflineFactory::buildOfflineObject(gridConfig, refinementConfig,
        regularizationConfig, densityEstimationConfig);
    offline->buildMatrix(grid.get(), regularizationConfig);
    offline->decomposeMatrix(regularizationConfig, densityEstimationConfig);
  }
  online = std::unique_ptr<DBMatOnlineDE>{DBMatOnlineDEFactory::buildDBMatOnlineDE(*offline,
     *grid, regularizationConfig.lambda_)};

  online->computeDensityFunction(alpha, newDataset, *grid,
      this->config->getDensityEstimationConfig(), true,
      this->config->getCrossvalidationConfig().enable_);
  online->setBeta(this->config->getLearnerConfig().beta);
  online->normalize(alpha, *grid);
  /*std::cout << getSurpluses().size() << std::endl;
  std::cout << getGrid().getSize() << std::endl;
  std::cout << getGrid().serialize() << std::endl;
  std::cout << getGrid().getDimension() << std::endl;*/
  // sgpp::datadriven::ModelFittingClassification test;
  // test.storeClassificator();

  // TEST
  std::string test = "[1,2,3,4]";
  test = "[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.011764705882352941, 0.01568627450980392, 0.023529411764705882, 0.047058823529411764, 0.07450980392156863, 0.09411764705882353, 0.10196078431372549, 0.09803921568627451, 0.07058823529411765, 0.043137254901960784, 0.03137254901960784, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01568627450980392, 0.027450980392156862, 0.10980392156862745, 0.3137254901960784, 0.47058823529411764, 0.5647058823529412, 0.592156862745098, 0.5725490196078431, 0.4745098039215686, 0.30196078431372547, 0.11764705882352941, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03137254901960784, 0.13333333333333333, 0.6470588235294118, 0.8352941176470589, 0.9294117647058824, 0.9607843137254902, 0.9686274509803922, 0.9607843137254902, 0.9333333333333333, 0.792156862745098, 0.3764705882352941, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1411764705882353, 0.5215686274509804, 1.0, 1.0, 1.0, 0.9882352941176471, 0.984313725490196, 0.9921568627450981, 1.0, 1.0, 0.9019607843137255, 0.35294117647058826, 0.10196078431372549, 0.0392156862745098, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3176470588235294, 0.788235294117647, 1.0, 0.996078431372549, 0.9254901960784314, 0.7725490196078432, 0.7450980392156863, 0.803921568627451, 0.9607843137254902, 1.0, 1.0, 0.792156862745098, 0.29411764705882354, 0.06666666666666667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.49019607843137253, 0.9294117647058824, 1.0, 0.9294117647058824, 0.5803921568627451, 0.20784313725490197, 0.19215686274509805, 0.23921568627450981, 0.6039215686274509, 0.9725490196078431, 1.0, 0.9294117647058824, 0.4823529411764706, 0.09411764705882353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.611764705882353, 0.9764705882352941, 1.0, 0.7607843137254902, 0.26666666666666666, 0.0, 0.0, 0.0, 0.19607843137254902, 0.7803921568627451, 0.9921568627450981, 0.9764705882352941, 0.6549019607843137, 0.1607843137254902, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6392156862745098, 0.9803921568627451, 1.0, 0.7411764705882353, 0.25098039215686274, 0.0, 0.0, 0.0, 0.11372549019607843, 0.611764705882353, 0.9686274509803922, 0.9921568627450981, 0.792156862745098, 0.24313725490196078, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5490196078431373, 0.9568627450980393, 1.0, 0.8470588235294118, 0.3607843137254902, 0.0, 0.0, 0.0, 0.06274509803921569, 0.4392156862745098, 0.9215686274509803, 1.0, 0.9254901960784314, 0.396078431372549, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3568627450980392, 0.7843137254901961, 0.9411764705882353, 0.7294117647058823, 0.28627450980392155, 0.0, 0.0, 0.0, 0.043137254901960784, 0.33725490196078434, 0.8705882352941177, 1.0, 0.9725490196078431, 0.596078431372549, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1411764705882353, 0.3568627450980392, 0.5098039215686274, 0.33725490196078434, 0.13333333333333333, 0.0, 0.0, 0.0, 0.03529411764705882, 0.25098039215686274, 0.8117647058823529, 1.0, 1.0, 0.8784313725490196, 0.011764705882352941, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.027450980392156862, 0.17254901960784313, 0.7490196078431373, 1.0, 1.0, 1.0, 0.10588235294117647, 0.00392156862745098, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.027450980392156862, 0.15294117647058825, 0.7176470588235294, 1.0, 1.0, 1.0, 0.2784313725490196, 0.00784313725490196, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023529411764705882, 0.1450980392156863, 0.7058823529411765, 1.0, 1.0, 1.0, 0.3058823529411765, 0.011764705882352941, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.027450980392156862, 0.1568627450980392, 0.7254901960784313, 1.0, 1.0, 1.0, 0.26666666666666666, 0.00784313725490196, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.027450980392156862, 0.1843137254901961, 0.7568627450980392, 1.0, 1.0, 1.0, 0.12549019607843137, 0.00392156862745098, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03137254901960784, 0.2235294117647059, 0.796078431372549, 1.0, 1.0, 0.9254901960784314, 0.027450980392156862, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0392156862745098, 0.054901960784313725, 0.08627450980392157, 0.09411764705882353, 0.08627450980392157, 0.043137254901960784, 0.023529411764705882, 0.0196078431372549, 0.043137254901960784, 0.30196078431372547, 0.8509803921568627, 1.0, 0.9882352941176471, 0.6901960784313725, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.10196078431372549, 0.27058823529411763, 0.4549019607843137, 0.48627450980392156, 0.43529411764705883, 0.24313725490196078, 0.08627450980392157, 0.03137254901960784, 0.058823529411764705, 0.41568627450980394, 0.9098039215686274, 1.0, 0.9294117647058824, 0.45098039215686275, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.30196078431372547, 0.7529411764705882, 0.9568627450980393, 0.9647058823529412, 0.9490196078431372, 0.8352941176470589, 0.3568627450980392, 0.07450980392156863, 0.09019607843137255, 0.5529411764705883, 0.9568627450980393, 0.996078431372549, 0.8509803921568627, 0.2901960784313726, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03529411764705882, 0.0784313725490196, 0.27058823529411763, 0.7686274509803922, 0.9803921568627451, 1.0, 1.0, 1.0, 0.996078431372549, 0.8196078431372549, 0.45098039215686275, 0.3254901960784314, 0.788235294117647, 0.9921568627450981, 0.9803921568627451, 0.6980392156862745, 0.17647058823529413, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.043137254901960784, 0.20392156862745098, 0.6980392156862745, 0.9803921568627451, 1.0, 1.0, 0.9372549019607843, 0.9764705882352941, 1.0, 0.9803921568627451, 0.8666666666666667, 0.7568627450980392, 0.9568627450980393, 1.0, 0.9333333333333333, 0.5058823529411764, 0.11764705882352941, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.058823529411764705, 0.36470588235294116, 0.9294117647058824, 1.0, 1.0, 1.0, 0.1803921568627451, 0.6588235294117647, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.796078431372549, 0.21176470588235294, 0.03529411764705882, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.07450980392156863, 0.43529411764705883, 0.9529411764705882, 1.0, 0.9450980392156862, 0.09411764705882353, 0.16862745098039217, 0.4549019607843137, 0.984313725490196, 1.0, 1.0, 1.0, 1.0, 1.0, 0.7725490196078432, 0.3411764705882353, 0.20784313725490197, 0.16470588235294117, 0.14901960784313725, 0.14901960784313725, 0.12156862745098039, 0.07058823529411765, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0784313725490196, 0.4627450980392157, 0.9568627450980393, 1.0, 0.9490196078431372, 0.48627450980392156, 0.596078431372549, 0.7764705882352941, 0.996078431372549, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9490196078431372, 0.8117647058823529, 0.7294117647058823, 0.6705882352941176, 0.6431372549019608, 0.6470588235294118, 0.5294117647058824, 0.23921568627450981, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06666666666666667, 0.39215686274509803, 0.9294117647058824, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.996078431372549, 0.9803921568627451, 0.9568627450980393, 0.9803921568627451, 0.996078431372549, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9372549019607843, 0.5333333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.043137254901960784, 0.2235294117647059, 0.6941176470588235, 0.9607843137254902, 1.0, 1.0, 1.0, 0.996078431372549, 0.9529411764705882, 0.8117647058823529, 0.5882352941176471, 0.3803921568627451, 0.615686274509804, 0.8117647058823529, 0.9137254901960784, 0.9647058823529412, 0.9921568627450981, 0.996078431372549, 1.0, 1.0, 0.9294117647058824, 0.5411764705882353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03529411764705882, 0.0784313725490196, 0.24705882352941178, 0.6470588235294118, 0.9137254901960784, 1.0, 0.9686274509803922, 0.8156862745098039, 0.5490196078431373, 0.33725490196078434, 0.17254901960784313, 0.07450980392156863, 0.1803921568627451, 0.3058823529411765, 0.4549019607843137, 0.596078431372549, 0.7137254901960784, 0.796078431372549, 0.8313725490196079, 0.796078431372549, 0.6549019607843137, 0.3254901960784314, 0.0, 0.0, 0.0]";
  DataVector dataVector = DataVector::fromString(test);
  std::cout << "works" << std::endl;
  std::cout << dataVector.toString() << std::endl;
  std::cout << "Evaluating" << std::endl;
  std::cout.precision(20);
  std::cout << "Prediction in MFDE: " << std::fixed << evaluate(dataVector) << std::endl;
}


bool ModelFittingDensityEstimationOnOff::refine(size_t newNoPoints,
    std::list<size_t> *deletedGridPoints) {
  // Coarsening, remove idx from alpha
  if (deletedGridPoints != nullptr && deletedGridPoints->size() > 0) {
    // Restructure alpha
    std::vector<size_t> idxToDelete{std::begin(*deletedGridPoints), std::end(*deletedGridPoints)};
    alpha.remove(idxToDelete);
  }
  // oldNoPoint refers to the grid size after coarsening
  auto oldNoPoints = alpha.size();

  // Refinement, expand alpha
  if (newNoPoints > oldNoPoints) {
    alpha.resizeZero(newNoPoints);
  }

  // Update online object: lhs, rhs and recompute the density function based on the b stored
  online->updateSystemMatrixDecomposition(config->getDensityEstimationConfig(),
      *grid, newNoPoints - oldNoPoints, *deletedGridPoints,
      config->getRegularizationConfig().lambda_);
  online->updateRhs(newNoPoints, deletedGridPoints);
  return true;
}

void ModelFittingDensityEstimationOnOff::update(Dataset& newDataset) {
  dataset = &newDataset;
  update(newDataset.getData());
}

void ModelFittingDensityEstimationOnOff::update(DataMatrix& newDataset) {
  if (grid == nullptr) {
    // Initial fitting of dataset
    fit(newDataset);
  } else {
    // Update the fit (streaming)
    online->computeDensityFunction(alpha, newDataset, *grid,
        this->config->getDensityEstimationConfig(), true,
        this->config->getCrossvalidationConfig().enable_);
    online->normalize(alpha, *grid);
  }
}

bool ModelFittingDensityEstimationOnOff::isRefinable() {
  if (grid != nullptr) {
    return online->getOfflineObject().isRefineable();
  }
  return false;
}

void ModelFittingDensityEstimationOnOff::reset() {
  grid.reset();
  online.reset();
  refinementsPerformed = 0;
}

}  // namespace datadriven
}  // namespace sgpp
