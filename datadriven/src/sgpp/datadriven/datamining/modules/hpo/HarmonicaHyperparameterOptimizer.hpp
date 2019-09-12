// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>

#include <memory>


namespace sgpp {
namespace datadriven {

/**
 * HarmonicaHyperparameterOptimizer coordinates data input, fitting and validation modules similarly to SparseGridMiner.
 * It offers access to different hyperparameter optimization procedures.
 */
class HarmonicaHyperparameterOptimizer : public HyperparameterOptimizer{
 public:
  /**
   * Constructor
   * @param miner configured instance of SGMiner object, that will provide the learning process.
   * The HyperparameterOptimizer instance will take ownership of the passed object.
   * @param fitterFactory configured instance of factory object that provides fitters with
   * manipulated hyperparameters. The HarmonicaHyperparameterOptimizer instance will take ownership of the passed object.
   * @param parser reference to parser object to read configuration info
   */
  HarmonicaHyperparameterOptimizer(SparseGridMiner* miner,
                          FitterFactory *fitterFactory,
                          DataMiningConfigParser &parser);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  HarmonicaHyperparameterOptimizer(const HarmonicaHyperparameterOptimizer &rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  HarmonicaHyperparameterOptimizer(HarmonicaHyperparameterOptimizer &&rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  HarmonicaHyperparameterOptimizer &operator=(HarmonicaHyperparameterOptimizer &&rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  HarmonicaHyperparameterOptimizer &operator=(const HarmonicaHyperparameterOptimizer &rhs) = delete;

  /**
   * Default destructor.
   */
  ~HarmonicaHyperparameterOptimizer() override = default;


  /**
   * Run hyperparameter optimization using Harmonica.
   */
  double run(bool writeToFile) override;
};
} /* namespace datadriven */
} /* namespace sgpp */
