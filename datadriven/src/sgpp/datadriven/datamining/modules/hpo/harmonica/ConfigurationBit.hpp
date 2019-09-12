// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ConfigurationBit_HPP
#define ConfigurationBit_HPP

#include <sgpp/datadriven/datamining/modules/hpo/harmonica/ConfigurationRestriction.hpp>

#include <string>
#include <vector>
#include <iostream>

namespace sgpp {
namespace datadriven {

/**
 * This class implements a boolean representation of hyperparameters for harmonica
 */
class ConfigurationBit {
 public:
  /**
   * Constructor
   * @param name to indentify the hyperparameter represented by this bit
   */
  explicit ConfigurationBit(std::string name)
      : name(name), constraints(0), value(0) {
    // std::cout << "Constructor: " << getName() << std::endl;
    // name = "test";
  }

  /**
   * Adds a reference to a new constraint limiting this bit
   * @param constraint
   */
  void addConstraint(ConfigurationRestriction *constraint);

  /**
   * removes last constraint (in case it was not valid)
   */
  void removeLastConstraint();

  /**
   * reset the value of this bit
   */
  void reset();

  /**
   * Set the value of this bit
   * @param input new value
   */
  void setValue(int input);

  /**
   * Get the value of this bit
   * @return
   */
  int getValue();

  /**
   * Get the name used to identify which hyperparameter is represented by this bit
   * @return name of the bit
   */
  std::string getName();

  void findComplexinner(std::string id, int bias);


 private:
  /**
   * name used to identify which hyperparameter is represented by this bit
   */
  std::string name;
  /**
   * vector pointing to the constraints that restrict this bit
   */
  std::vector<ConfigurationRestriction *> constraints;
  /**
   * current value of this bit (-1 or 1 or 0 for unset)
   */
  int value;
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* ConfigurationBit_HPP */
