// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_Harmonica_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_Harmonica_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>

#include <vector>
#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Class to host all methods required to perform the harmonica algorithm
 */
class Harmonica {
 public:
  /**
   * Constructor
   * @param fitterFactory to produce fitter type objects
   */
  explicit Harmonica(FitterFactory *fitterFactory);

  /**
   * First step in harmonica. Configurations are prepared for evaluation and the parity
   * function matrix is constructed for later use.
   * @param fitters container to store fitters for evaluation outside the class
   * @param seed for random sampling
   * @param configStrings container to store information about the configurations in string form
   */
  void prepareConfigs(std::vector<ModelFittingBase*> &fitters,
                      int seed,
                      std::vector<std::string> &configStrings);

  /**
   * Function to create a vector of random numbers within the valid range of possible configurations
   * without duplicates
   * @param nBits number of bits a configuration consists of
   * @param configIDs container to store generated random configuration ID's
   * @param seed to use in random generator
   * @param start to offset newly generated ID's from existing ones from previous iterations
   */
  void createRandomConfigs(size_t nBits, std::vector<int> &configIDs, int seed, size_t start);
  /**
   * Second step of the harmonica algorithm. Calculates relavance of ConfigurationBits and introduces
   * constraints to reduce the search space
   * @param transformedScores input (possibly transformed) scores after evaluation
   * @param lambda used for regression
   * @param shrink number of constraints to introduce
   */
  void calculateConstrainedSpace(const DataVector &transformedScores, double lambda, int shrink);
  /**
   * Transforms scores to accentuate the optimum
   * @param source
   * @param target
   */
  void transformScores(const DataVector &source, DataVector &target);
  /**
   * resolves constraints, fixing free and dependent bits
   */
  void fixConfigBits(bool resetFree);
  /**
   * resets bits to be able to resolve constraints again
   */
  void resetBits();
  /**
   * Sets bits and the resulting parameter configuration based on a configID while simultaneously
   * creating one row of the parity function matrix
   * @param configID bits of the configuration as an integer
   * @param matrixrow index of the row of the parity matrix to be filled
   */
  void setParameters(int configID, size_t matrixrow);
  /**
   * Adds a constraint based on an entry in the parity row.
   * @param idx index of the parity row that holds pointers to the bits that will be constrained
   * @param bias bias value of the constraint
   */
  bool addConstraint(size_t idx, int bias);
  /**
   * Tests all constraints for validity in the current bit configuration
   * @return true when all constraints are met
   */
  bool checkConstraints();
  /**
   * Recalculates a configID to the new constrained binary space
   * @param configID binary configuration in old space as integer
   * @param oldFreeBits bits that received the old binary values
   * @return binary configuartion in new space as integer
   */
  int moveToNewSpace(int configID, std::vector<ConfigurationBit *> oldFreeBits);

 protected:
  /**
   * matrix that holds the values of the parity function
   * (predictors for regression) for all samples
   */
  base::DataMatrix paritymatrix;
  /**
   * pointer to the fitterFactory to produce fitters for each configuration
   */
  FitterFactory *fitterFactory;
  /**
   * binary configurations as integers
   */
  std::vector<int> configIDs;
  /**
   * Scores saved for moving configurations to new space
   */
  DataVector savedScores;
  /**
   * row to create the parity matrix and the constraints
   */
  std::vector<std::vector<ConfigurationBit *> > parityrow;
  /**
   * bits that are not dependent on other bits
   */
  std::vector<ConfigurationBit *> freeBits;
  /**
   * all configuration bits
   */
  std::vector<ConfigurationBit *> configBits;
  /**
   * all constraints that currently exist
   */
  std::vector<std::unique_ptr<ConfigurationRestriction>> constraints;
};
} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_Harmonica_HPP_ */
