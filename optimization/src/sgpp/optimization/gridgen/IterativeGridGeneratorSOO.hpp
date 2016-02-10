// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP

#include <cstddef>
#include <functional>

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>

namespace SGPP {
namespace optimization {

/**
 * Iterative grid generator using simultaneous optimistic
 * optimization (SOO).
 */
class IterativeGridGeneratorSOO : public IterativeGridGenerator {
 public:
  typedef std::function<size_t(size_t)> AdaptivityFunction;

  /// default adaptivity
  static constexpr float_t DEFAULT_ADAPTIVITY = 0.5;

  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param f             objective function
   * @param grid          grid (should be empty)
   * @param N             maximal number of grid points
   * @param adaptivity    adaptivity (positive number)
   */
  IterativeGridGeneratorSOO(ScalarFunction& f,
                            base::Grid& grid,
                            size_t N,
                            float_t adaptivity = DEFAULT_ADAPTIVITY);

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
  void setAdaptivity(float_t adaptivity);

  /*
   * @param adaptivity  adaptivity (function of the form
   *                    "iteration number --> maximal refinement depth")
   */
  void setAdaptivity(AdaptivityFunction adaptivity);

 protected:
  /// adaptivity
  AdaptivityFunction hMax;
};

}
}

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORSOO_HPP */
