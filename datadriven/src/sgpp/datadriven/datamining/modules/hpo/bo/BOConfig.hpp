// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CLION_BOCONFIG_HPP
#define CLION_BOCONFIG_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <vector>
#include <random>

namespace sgpp {
namespace datadriven {

/**
 * Container class to store a conrete hyperparameter configuration for interaction with Bayesian Optimization
 */
class BOConfig {
 public:
  /**
   * Default Constructor.
   */
  BOConfig() = default;

  /**
   * Constructor for making a prototype based on the number of hyperparameters
   * @param discOptions number of options for each discrete parameter
   * @param catOptions number of options for each categorical parameter
   * @param nCont number of continuous parameters
   */
  BOConfig(std::vector<int> *discOptions, std::vector<int> *catOptions, size_t nCont);

  /**
   * Iterator over discrete parameter options
   * @return stopping criterion
   */
  bool nextDisc();

  /**
   * calculation of discrete part of the distance between two BOConfigs/sample points
   * @param other sample point to calculate distance to
   * @param scales scaling of hyperparameters in relation to each other
   */
  void calcDiscDistance(BOConfig &other, base::DataVector &scales);

  /**
   * finish previous distance calculation by adding the continuous part
   * @param input continuous part of the other (new) sample point
   * @param scales scaling of hyperparameters in relation to each other
   * @return complete distance measure
   */
  double getTotalDistance(const base::DataVector &input, base::DataVector &scales);

  /**
   * Get number of continuous parameters
   * @return number of continuous parameters
   */
  size_t getContSize();

  /**
 * Get number total number of parameters
 * @return number of continuous parameters
 */
  size_t getNPar() const;

  /**
   * Set the continuous parameters according to input
   * @param input DataVector holding continuous parameters
   */
  void setCont(const base::DataVector &input);

  /**
   * Get the value of a specific continuous parameter
   * @param idx parameter position
   * @return parameter value
   */
  double getCont(size_t idx);
  /**
   * Get the value of a specific discrete parameter
   * @param idx parameter position
   * @return parameter value
   */
  int getDisc(size_t idx);
  /**
   * Get the value of a specific categorical parameter
   * @param idx parameter position
   * @return parameter value
   */
  int getCat(size_t idx);

  /**
   * Set score measured on this sample
   * @param input score
   */
  void setScore(double input);

  /**
   * Get score measured on this sample
   * @return score
   */
  double getScore();

  /**
   * Compute complete distance to another BOConfig/sample point
   * @param other sample point to calculate distance to
   * @param scales scaling of hyperparameters in relation to each other
   * @return distance measure
   */
  double getScaledDistance(BOConfig &other, const base::DataVector &scales);

  /**
   * Generate a random config
   * @param generator for seeded rng
   */
  void randomize(std::mt19937 &generator);

 private:
  /**
   * DataVector containing values of continuous parameters (shifted to [0,1])
   */
  base::DataVector cont;
  /**
   * Vector containing values of discrete parameters
   */
  std::vector<int> disc;
  /**
   * Vector containing values of categorical parameters
   */
  std::vector<int> cat;
  /**
   * pointer to vector containing number of options for each discrete parameter
   */
  std::vector<int> *discOptions = nullptr;
  /**
   * pointer to vector containing number of options for each categorical parameter
   */
  std::vector<int> *catOptions = nullptr;
  /**
   * score measured when evaluating this config
   */
  double score = 0;
  /**
   * discrete component of the distance to some potential new sample point
   */
  double discDistance = 0;
};
} /* namespace datadriven */
} /* namespace sgpp */

#endif   // CLION_BOCONFIG_HPP
