/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMiner.hpp
 *
 * Created on: Oct 7, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/Visualizer.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <initializer_list>
#include <memory>
#include <sstream>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * SparseGridMiner models the entire mining process for data mining with sparse grids. It aggregates
 * and automates data input, fitting and validation modules and controls the mining process.
 */
class SparseGridMiner {
 public:
  /**
   * Constructor
   * @param fitter configured instance of fitter object that generalize the model. The miner
   * instance will take ownership of the passed object.
   * @param scorer configured instance of scorer object that will assess the quality of the
   * generalization provided by the fitter on testing data. The miner instance will take ownership
   * of the passed object.
   * @param visualizer configured instance of the visualizer object that will execute the
   * visualization module of the model.
   * The miner instance will take ownership of the passed object
   */
  SparseGridMiner(ModelFittingBase *fitter, Scorer *scorer, Visualizer *visualizer);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  SparseGridMiner(const SparseGridMiner &rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  SparseGridMiner(SparseGridMiner &&rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  SparseGridMiner &operator=(SparseGridMiner &&rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  SparseGridMiner &operator=(const SparseGridMiner &rhs) = delete;

  /**
   * Default destructor.
   */
  virtual ~SparseGridMiner() = default;

  /**
   * Perform Learning cycle: Get samples from data source and based on the scoring procedure,
   * generalize data by fitting and asses quality of the fit.
   */
  virtual double learn(bool verbose) = 0;

  /**
   * Returns the trained model
   * @return the trained model
   */
  ModelFittingBase *getModel();

  void setModel(ModelFittingBase *model);

  Visualizer *getVisualizer();

  /**
   * Evaluate the model on a certain test dataset.
   *
   * @param testDataset dataset used quantify accuracy using #sgpp::datadriven::Metric.
   * @return score of the fit.
   */

  double test(Dataset &testDataset);

  /**
   * Print output on one process.
   * @param message
   */
  static void print(const std::string &message);

  /**
   * Print output on one process.
   * @param message
   */
  static void print(const char *message);

  /**
   * Print output on one process.
   * @param messageStream stream with the concatenated message
   */
  static void print(std::ostringstream &messageStream);

 protected:
  /**
   * Fitter that trains a model based on data samples.
   */
  std::unique_ptr<ModelFittingBase> fitter;
  /**
   * Scorer that quantifies the quality of a fit. (e.g. cross validation or training with testing)
   */
  std::unique_ptr<Scorer> scorer;

 /*
  * Visualizer which generates files as
  * an input of graphic libraries to visualize the models
  */
  std::unique_ptr<Visualizer> visualizer;
};
}  // namespace datadriven
}  // namespace sgpp
