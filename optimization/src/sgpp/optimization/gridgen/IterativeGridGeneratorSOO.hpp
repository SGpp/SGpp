// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

#include <cstddef>
#include <functional>

namespace sgpp {
namespace optimization {

/**
 * Iterative grid generator using simultaneous optimistic
 * optimization (SOO).
 */
class IterativeGridGeneratorSOO : public IterativeGridGenerator {
 public:
  typedef std::function<size_t(size_t)> AdaptivityFunction;

  /// default adaptivity
  static constexpr double DEFAULT_ADAPTIVITY = 0.5;

  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param f             objective function
   * @param grid          grid (should be empty)
   * @param N             maximal number of grid points
   * @param adaptivity    adaptivity (positive number)
   */
  IterativeGridGeneratorSOO(base::ScalarFunction& f, base::Grid& grid, size_t N,
                            double adaptivity = DEFAULT_ADAPTIVITY);

  /**
   * Destructor.
   */
  ~IterativeGridGeneratorSOO() override;

  /**
   * Generate the grid.
   *
   * @return true on success, otherwise false
   */
  bool generate() override;

  /*
   * @return            adaptivity (function of the form
   *                    "iteration number --> maximal refinement depth")
   */
  AdaptivityFunction getAdaptivity() const;

  /*
   * @param adaptivity  adaptivity (positive number)
   */
  void setAdaptivity(double adaptivity);

  /*
   * @param adaptivity  adaptivity (function of the form
   *                    "iteration number --> maximal refinement depth")
   */
  void setAdaptivity(AdaptivityFunction adaptivity);

 protected:
  /// adaptivity
  AdaptivityFunction hMax;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP */
