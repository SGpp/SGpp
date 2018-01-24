// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsebase::Grids.org

#ifndef ConfigurationBit_HPP
#define ConfigurationBit_HPP

#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationRestriction.hpp>

#include <sgpp/globaldef.hpp>
#include <list>

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
  ConfigurationBit()
      : constraints(), value(0), bVisited(false) {}

  /**
   * Destructor
   */
  // ~ConfigurationBit() {}

  void addConstraint(ConfigurationRestriction* constraint);
  
  int evaluate(int* input);
  
  void reset();
  // void mult(base::DataVector& alpha, base::DataVector& result);
  // void multTranspose(base::DataVector& source, base::DataVector& result);

  // double getDuration();

 protected:
  /// reference to the base::Grid's base::GridStorage object
  std::list<ConfigurationRestriction> constraints;
  int value;
  bool bVisited;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* ConfigurationBit_HPP */
