// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsebase::Grids.org

#ifndef ConfigurationRestriction_HPP
#define ConfigurationRestriction_HPP

// #include <sgpp/datadriven/datamining/modules/hpo/ConfigurationBit.hpp>


// EDIT: cyclic dependency


#include <sgpp/globaldef.hpp>
#include <vector>

namespace sgpp {
namespace datadriven {
class ConfigurationBit;

/**
 * This class implements OperationB for a base::Grids with linear basis ansatzfunctions without boundaries
 */
class ConfigurationRestriction {
 public:
  /**
   * Constructor of OperationBLinear
   *
   * @param base::Grid base::Grid
   * @param dataset the dataset that should be evaluated
   */
  ConfigurationRestriction(std::vector<ConfigurationBit*> &parameters, int bias);

  /**
   * Destructor
   */
  // ~ConfigurationRestriction() {}



  void reduceOpenBits();

  int getOpenBits();

  void resolve();

  bool check();

  void reset();

 protected:
  /// reference to the base::Grid's base::GridStorage object
  std::vector<ConfigurationBit*> parameters;
  int bias;
  int openBits = 0;
  // bool visited;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* ConfigurationRestriction_HPP */
