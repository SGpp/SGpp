// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsebase::Grids.org

#ifndef ConfigurationBit_HPP
#define ConfigurationBit_HPP

#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationRestriction.hpp>

#include <sgpp/globaldef.hpp>
#include <list>
#include <vector>
#include <iostream>

namespace sgpp {
namespace datadriven {

/**
 * This class implements OperationB for a base::Grids with linear basis ansatzfunctions without boundaries
 */
class ConfigurationBit {
 public:
  /**
   * Constructor of OperationBLinear
   *
   * @param base::Grid base::Grid
   * @param dataset the dataset that should be evaluated
   */
  explicit ConfigurationBit(std::string namen)
      :name(namen), constraints(0), value(0){
    //std::cout << "Constructor: " << getName() << std::endl;
    //name = "test";
  }

  //~ConfigurationBit(){
  //  std::cout << "ConfigurationBit destroyed!" << name<<std::endl;
  //}

  /**
   * Destructor
   */
  // ~ConfigurationBit() {}

  void addConstraint(ConfigurationRestriction* constraint);

  void removeLastConstraint();

  void reset();

  void setValue(int input);

  int getValue();

  std::string getName();
  // void mult(base::DataVector& alpha, base::DataVector& result);
  // void multTranspose(base::DataVector& source, base::DataVector& result);

  // double getDuration();


 private:
  std::string name;
  std::vector<ConfigurationRestriction*> constraints;
  int value;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* ConfigurationBit_HPP */
