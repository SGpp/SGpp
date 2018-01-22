// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsebase::Grids.org

#ifndef ConfigurationRestriction_HPP
#define ConfigurationRestriction_HPP


#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace datadriven {

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
  ConfigurationRestriction(std::list<int> parameters, int bias)
      : parameters(parameters), bias(bias) {}

  /**
   * Destructor
   */
  // ~ConfigurationRestriction() {}

  // void mult(base::DataVector& alpha, base::DataVector& result);
  // void multTranspose(base::DataVector& source, base::DataVector& result);

  // double getDuration();

 protected:
  /// reference to the base::Grid's base::GridStorage object
  std::list<int> parameters;
  int bias;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* ConfigurationRestriction_HPP */
