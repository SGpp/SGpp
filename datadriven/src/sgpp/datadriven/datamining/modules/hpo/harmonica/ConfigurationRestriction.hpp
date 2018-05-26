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
 * This class represents a constraint on the ConfigurationBits used by harmonica
 */
class ConfigurationRestriction {
 public:
  /**
   * Constructor
   * @param parameters bits restricted by this constraint
   * @param bias value the product of the restricted bits must be equal to
   */
  ConfigurationRestriction(std::vector<ConfigurationBit*> &parameters, int bias);

  /**
   * Decrease the counter of open bits during constraint resolution
   */
  void reduceOpenBits();

  /**
   * get counter of open bits during constraint resolution
   * @return open bit counter
   */
  int getOpenBits();

  /**
   * resolve constraint when there is only one open bit remaining
   */
  void resolve();

  /**
   * checks whether this constraint is satisfied with the current bit configuration
   * @return true when constraint is satisfied
   */
  bool check();

  /**
   * resets constraint resolution state
   */
  void reset();

 protected:
  /**
   * bits affected by the constraint
   */
  std::vector<ConfigurationBit*> parameters;
  /**
   * value the product of the restricted bits must be equal to
   */
  int bias;
  /**
   * number of open (unset) bits during constraint resolution
   */
  int openBits = 0;
  // bool visited;
};

}  // namespace datadriven
}  // namespace sgpp

#endif /* ConfigurationRestriction_HPP */
