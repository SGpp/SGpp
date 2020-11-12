// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/exception/not_implemented_exception.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BOConfig.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/bo/BayesianOptimization.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>
#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/harmonica/Harmonica.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/BoHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HarmonicaHyperparameterOptimizer.hpp>
#include <sgpp/datadriven/datamining/builder/LeastSquaresRegressionMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationLeastSquares.hpp>

#include <string>
#include <vector>

using sgpp::datadriven::Dataset;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::datadriven::BOConfig;
using sgpp::datadriven::ConfigurationBit;

BOOST_AUTO_TEST_SUITE(HPOTest)

class ModelFittingTester : public sgpp::datadriven::ModelFittingBase {
 public:
  ModelFittingTester(double x, int y, int z) {
    value = x * x + y * y / 100.0 - z / 1.4;
    if (z == 2) {
      value *= -1;
    }
    value += 1;
    if (value < 0) {
      std::cout << "Error! value < 0" << std::endl;
    }
    config.reset(new sgpp::datadriven::FitterConfigurationLeastSquares());
  }

  void fit(Dataset &dataset) override {}

  void fit(Dataset &, Dataset &) override {}

  bool adapt() override { return false; }

  void update(Dataset &dataset) override {}

  void update(Dataset &, Dataset &) override {}

  void reset() override {}

  void resetTraining() override {}

  double evaluate(const DataVector &sample) override { return -420; }

  void evaluate(DataMatrix &samples, DataVector &results) override { results[0] = sqrt(value); }

  double computeResidual(DataMatrix &validationData) const override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingTester::computeResidual() not implemented!");
  }

  void updateRegularization(double lambda) override {
    throw sgpp::base::not_implemented_exception(
        "ModelFittingTester::updateRegularization() not implemented!");
  }

  double value;
};

class FitterFactoryTester : public sgpp::datadriven::FitterFactory {
 public:
  FitterFactoryTester() {
    // create parameters
    conpar["x"] = sgpp::datadriven::ContinuousParameter(4, "x", -1, 1);
    dispar["y"] = sgpp::datadriven::DiscreteParameter("y", -10, 10);
    catpar["c"] = sgpp::datadriven::DiscreteParameter("c", 0, 2);
  }

  sgpp::datadriven::ModelFittingBase *buildFitter() override {
    // making model from parameter values directly
    return new ModelFittingTester(conpar["x"].getValue(), dispar["y"].getValue(),
                                  catpar["c"].getValue());
  }
};

class FitterFactoryTesterHarm : public sgpp::datadriven::FitterFactory {
 public:
  FitterFactoryTesterHarm() {
    for (int i = 0; i < 12; ++i) {
      exconfBits.emplace_back(std::to_string(i));
    }
  }

  sgpp::datadriven::ModelFittingBase *buildFitter() override { return nullptr; }

  void getConfigBits(std::vector<ConfigurationBit *> &configBits) override {
    for (auto &exconfBit : exconfBits) {
      configBits.push_back(&exconfBit);
    }
  }

  std::vector<ConfigurationBit> exconfBits;
};

class HarmonicaTester : public sgpp::datadriven::Harmonica {
 public:
  explicit HarmonicaTester(sgpp::datadriven::FitterFactory *fft) : Harmonica(fft) {}
  std::vector<std::vector<ConfigurationBit *>> &getParityrow() { return parityrow; }
  std::vector<ConfigurationBit *> getFreeBits() { return freeBits; }
};

/*
 * To fix: Verbose miner
 * dummy set account for validation
 */

BOOST_AUTO_TEST_CASE(upperLevelTest) {
  // using actual files for (dummy) data and config
  std::string path("datadriven/tests/pipeline/config_hpo.json");
  sgpp::datadriven::DataMiningConfigParser parser(path);
  // sgpp::datadriven::DataSourceBuilder builder;
  // sgpp::datadriven::DataSourceConfig config;
  // parser.getDataSourceConfig(config, config);
  sgpp::datadriven::LeastSquaresRegressionMinerFactory minfac{};

  sgpp::datadriven::BoHyperparameterOptimizer bohpo(minfac.buildMiner(path),
                                                    new FitterFactoryTester(), parser);
  sgpp::datadriven::HarmonicaHyperparameterOptimizer harmhpo(minfac.buildMiner(path),
                                                             new FitterFactoryTester(), parser);
  double res1 = bohpo.run(false);
  double res2 = harmhpo.run(false);
  // testing arbitrary performance lower bound
  // theoretical optimum is 0.4/1.4 ~= 0.2857 for c = 1, x = 0, y = 0
  BOOST_CHECK_LE(res1, 0.3);
  BOOST_CHECK_LE(res2, 0.3);
}

BOOST_AUTO_TEST_CASE(harmonicaConfigs) {
  // tests the bit management, especially setParameters and addConstraint by comparing to a vector
  // of all possible bit configurations
  std::mt19937 generator(34);
  FitterFactoryTesterHarm fft{};
  HarmonicaTester harmonica(&fft);
  bool testar[4096];
  int oldidar[4096];
  int nIDs = 4096;
  std::vector<sgpp::datadriven::ModelFittingBase *> fitters(1);
  std::vector<std::string> configStrings(1);
  std::vector<sgpp::datadriven::ConfigurationRestriction> constraints{};

  for (int i = 0; i < 3; ++i) {
    harmonica.prepareConfigs(fitters, 33, configStrings);
    std::uniform_int_distribution<int> distcons(
        0, static_cast<int>(harmonica.getParityrow().size() - 1));
    std::geometric_distribution<int> distgeo(0.05 + 0.2 * i);  // not completely safe but okay
    std::uniform_int_distribution<int> distbias(0, 1);
    for (int k = 0; k < 4096; ++k) {
      oldidar[k] = -1;
    }
    std::vector<ConfigurationBit *> freeBits = harmonica.getFreeBits();
    for (int j = 0; j < nIDs; ++j) {
      harmonica.setParameters(j, 0);
      int u = 0;
      int x = 1;
      for (auto &bit : fft.exconfBits) {
        u = u + x * ((bit.getValue() + 1) / 2);
        x = x * 2;
      }
      oldidar[u] = j;
    }
    for (int l = 0; l < 5; ++l) {
      int consid;
      if (l % 2 == 0) {
        consid = distcons(generator);
      } else {
        consid = distgeo(generator);
      }
      int bias = -1 + 2 * distbias(generator);

      constraints.emplace_back(harmonica.getParityrow()[consid], bias);
      /*
       std::cout << "New Bit: ";
       for(auto &bit : harmonica.getParityrow()[consid]){
       std::cout << bit->getName() << ",";
       }
       std::cout << "Bias: " << bias << std::endl;
       */
      int cnttrue = 0;

      for (int k = 0; k < 4096; ++k) {
        testar[k] = false;
      }
      for (int j = 0; j < nIDs; ++j) {
        harmonica.setParameters(j, 0);
        int v = 0;
        int m = 1;
        for (auto &bit : fft.exconfBits) {
          v = v + m * ((bit.getValue() + 1) / 2);
          m = m * 2;
        }
        BOOST_CHECK_EQUAL(harmonica.moveToNewSpace(oldidar[v], freeBits), j);
        BOOST_CHECK(!testar[v]);  // duplicate parameter set
        testar[v] = true;
        bool valid = true;
        for (auto &constraint : constraints) {
          if (!constraint.check()) {
            valid = false;
          }
        }
        if (valid) {
          cnttrue++;
        }
      }
      bool added = harmonica.addConstraint(consid, bias);  // if not added, freebits broken
                                                           /*
                                                            std::cout << "Counter true: " << cnttrue << " | Added: "<< added
                                                            <<std::endl;
                                                            for (auto &bit : fft.exconfBits) {
                                                            std::cout << bit.getName() << ": " << bit.getValue() << ", ";
                                                            }
                                                            std::cout << std::endl;
                                                            */
      // either all, none or half are viable
      BOOST_CHECK(((!added) == (cnttrue == 0)));
      BOOST_CHECK((added == (cnttrue == nIDs || 2 * cnttrue == nIDs)));
      if (added) {
        nIDs = cnttrue;
      } else {
        harmonica.resetBits();
        harmonica.fixConfigBits(true);
        constraints.pop_back();
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(DistanceCalculation) {
  // Comparing the two methods of calculating distance in bayesian optimization.
  // Should be equal
  std::mt19937 generator(123);
  std::uniform_real_distribution<double> rand(0.0, 1.0);
  std::vector<int> discOptions = {2, 3};
  std::vector<int> catOptions = {2, 3};
  size_t nCont = 7;
  BOConfig prototype{&discOptions, &catOptions, nCont};
  DataVector scales(prototype.getNPar() + 1, 1);

  BOConfig oldconf(prototype);
  oldconf.randomize(generator);
  for (int j = 0; j < 100; ++j) {
    for (size_t i = 0; i < scales.size(); ++i) {
      scales[i] = rand(generator);
    }
    BOConfig curconf(prototype);
    curconf.randomize(generator);
    oldconf.calcDiscDistance(curconf, scales);
    double curdist = oldconf.getScaledDistance(curconf, scales);
    DataVector contpoint(nCont);
    for (size_t k = 0; k < nCont; ++k) {
      contpoint[k] = curconf.getCont(k);
    }
    double twostepdist = oldconf.getTotalDistance(contpoint, scales);
    BOOST_CHECK_CLOSE(curdist, twostepdist, 1e-5);
    oldconf = curconf;
  }
}

BOOST_AUTO_TEST_CASE(addSamplesGP) {
  // Adding samples to Gaussian Process and correctly reading out their value
  std::vector<BOConfig> initialConfigs{};
  std::mt19937 generator(123);

  std::vector<int> discOptions = {2, 3};
  std::vector<int> catOptions = {2, 3};
  size_t nCont = 2;
  BOConfig prototype{&discOptions, &catOptions, nCont};

  std::vector<double> scores = {0, 42, 21, 30, 5};
  initialConfigs.reserve(scores.size());

  for (size_t i = 0; i < scores.size(); i++) {
    initialConfigs.emplace_back(prototype);
    initialConfigs[i].randomize(generator);
    initialConfigs[i].setScore(scores[i]);
  }

  // single point not possible because of normalize
  sgpp::datadriven::BayesianOptimization bo(initialConfigs);

  DataVector scales(initialConfigs.front().getNPar() + 1, 1);

  DataMatrix kernelmatrix(initialConfigs.size(), initialConfigs.size());

  double noise = pow(10, -scales.back() * 10);
  for (size_t i = 0; i < initialConfigs.size(); ++i) {
    for (size_t k = 0; k < i; ++k) {
      double tmp = bo.kernel(initialConfigs[i].getScaledDistance(initialConfigs[k], scales));
      kernelmatrix.set(k, i, tmp);
      kernelmatrix.set(i, k, tmp);
    }
    kernelmatrix.set(i, i, 1 + noise);
  }

  DataVector kernelrow(initialConfigs.size());

  DataVector dscores(scores);
  dscores.normalize();
  dscores.sub(DataVector(dscores.size(), dscores.sum() / static_cast<double>(dscores.size())));

  for (size_t i = 0; i < initialConfigs.size(); i++) {
    kernelmatrix.getColumn(i, kernelrow);
    BOOST_CHECK_CLOSE(bo.mean(kernelrow), dscores[i], 1e-5);
    BOOST_CHECK_CLOSE(bo.var(kernelrow, 1), 0, 1e-5);
    // std::cout << "Mean " << i << ": " << bo.mean(kernelrow) << "  |
    // Original: "
    // << dscores[i] << " | Var: " << bo.var(kernelrow, 1) <<std::endl;
  }

  scores = {0, 42, 21, 30, 5, 33, 11, 15, 4.6};

  for (size_t i = 5; i < scores.size(); ++i) {
    initialConfigs.emplace_back(prototype);
    initialConfigs[i].randomize(generator);
    initialConfigs[i].setScore(scores[i]);
    bo.updateGP(initialConfigs[i], true);
  }

  kernelmatrix.resize(initialConfigs.size(), initialConfigs.size());
  kernelrow.resize(initialConfigs.size());
  dscores = DataVector(scores);

  for (size_t i = 0; i < initialConfigs.size(); ++i) {
    for (size_t k = 0; k < i; ++k) {
      double tmp = bo.kernel(initialConfigs[i].getScaledDistance(initialConfigs[k], scales));
      kernelmatrix.set(k, i, tmp);
      kernelmatrix.set(i, k, tmp);
    }
    kernelmatrix.set(i, i, 1 + noise);
  }

  dscores.normalize();
  dscores.sub(DataVector(dscores.size(), dscores.sum() / static_cast<double>(dscores.size())));

  for (size_t i = 0; i < initialConfigs.size(); i++) {
    kernelmatrix.getColumn(i, kernelrow);
    BOOST_CHECK_CLOSE(bo.mean(kernelrow), dscores[i], 1e-5);
    BOOST_CHECK_CLOSE(bo.var(kernelrow, 1), 0, 1e-5);
    // std::cout << "Mean " << i << ": " << bo.mean(kernelrow) << "  |
    // Original: " << dscores[i] << " | Var: " << bo.var(kernelrow, 1) <<std::endl;
  }
}

BOOST_AUTO_TEST_CASE(fitScalesGP) {
  // test gaussian process fitting by fitting to a second GP
  std::vector<BOConfig> initialConfigs{};
  std::mt19937 generator(100);

  std::vector<int> discOptions = {2, 3};
  std::vector<int> catOptions = {2, 3};
  size_t nCont = 6;
  BOConfig prototype{&discOptions, &catOptions, nCont};

  DataVector scores{std::vector<double>({-0.1, 0.1})};

  initialConfigs.reserve(scores.size());

  for (size_t i = 0; i < scores.size(); i++) {
    initialConfigs.emplace_back(prototype);
    initialConfigs[i].randomize(generator);
    initialConfigs[i].setScore(scores[i]);
  }

  // single point not possible because of normalize
  sgpp::datadriven::BayesianOptimization genprocess(initialConfigs);
  sgpp::datadriven::BayesianOptimization fitprocess(initialConfigs);

  DataVector scales{1, 0.1, 0.3, 0.4, 0.01, 0.9, 0.2, 0.8, 1, 0.5, 0.3};

  scales.resize(prototype.getNPar() + 1);
  genprocess.setScales(scales, 1);
  DataVector avscales(scales.size(), 1);

  for (int j = 0; j < 100; ++j) {
    DataVector kernelrow(initialConfigs.size());
    BOConfig npoint(prototype);
    npoint.randomize(generator);
    for (size_t i = 0; i < initialConfigs.size(); i++) {
      kernelrow[i] = genprocess.kernel(npoint.getScaledDistance(initialConfigs[i], scales));
    }
    double mean = genprocess.mean(kernelrow);
    double var = genprocess.var(kernelrow, 1.001);
    // std::cout << "Mean: " << mean << " | sqVar: " << sqrt(var) << std::endl;
    std::normal_distribution<double> distribution(mean, var);
    npoint.setScore(distribution(generator));
    initialConfigs.push_back(npoint);

    scores.push_back(npoint.getScore());
    genprocess.updateGP(npoint, false);
    fitprocess.updateGP(npoint, false);
    DataVector fitscales(fitprocess.fitScales());
    // std::cout << fitscales.toString() << std::endl;
    // fitscales.sub(scales);
    // std::cout << "Diff: " << fitscales.l2Norm()
    // << " | Max Norm: " << fitscales.maxNorm() << std::endl;
    // fitscales.add(scales);
    fitscales.mult(0.1);
    avscales.mult(0.9);
    avscales.add(fitscales);
    avscales.sub(scales);
    // std::cout << "################################ AvDiff: " <<
    // avscales.l2Norm()
    // << " | Max Norm: " << avscales.maxNorm() << " | size ratio: " <<
    // sizeratio << std::endl;
    if (j > 90) {
      // difference is allowed to rise with the squareroot of the dimensionality
      BOOST_CHECK_LE(avscales.l2Norm(), sqrt(static_cast<double>(prototype.getNPar() + 1)) * 0.25);
    }
    avscales.add(scales);
  }
}

BOOST_AUTO_TEST_CASE(validAcquisitionFunction) {
  // testing acquisition function for monotonicity with respect to mean and
  // variance
  // not every acquisition function fullfills this but expected improvement does
  std::mt19937 generator(123);
  std::uniform_real_distribution<double> ranvar(0.0, 1.0);
  std::uniform_real_distribution<double> ranmean(-1.0, 1.0);
  double oldvar = ranvar(generator);
  double oldmean = ranmean(generator);
  double oldbest = ranmean(generator);
  double oldAc = sgpp::datadriven::BayesianOptimization::acquisitionEI(oldmean, oldvar, oldbest);
  for (int j = 0; j < 1000; ++j) {
    double curvar = ranvar(generator);
    double curmean = ranmean(generator);
    double curbest = ranmean(generator);
    double curAc = sgpp::datadriven::BayesianOptimization::acquisitionEI(curmean, curvar, oldbest);
    if (curmean < oldmean && curvar > oldvar) {
      BOOST_CHECK_LT(curAc, oldAc);
      if (oldAc < curAc) {
        std::cout << curmean << "," << curvar << "," << curbest << "," << curAc << std::endl;
        std::cout << oldmean << "," << oldvar << "," << oldbest << "," << oldAc << std::endl;
        std::cout << std::endl;
      }
    }
    oldvar = curvar;
    oldmean = curmean;
    oldbest = curbest;
    oldAc = sgpp::datadriven::BayesianOptimization::acquisitionEI(curmean, curvar, curbest);
  }
}

BOOST_AUTO_TEST_SUITE_END()
