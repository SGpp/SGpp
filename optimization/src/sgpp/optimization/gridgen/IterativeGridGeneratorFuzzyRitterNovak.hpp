// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORFUZZYRITTERNOVAK_HPP
#define SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORFUZZYRITTERNOVAK_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

#include <cstddef>
#include <vector>

namespace sgpp {
namespace optimization {

/**
 * Iterative grid generation based on Ritter/Novak's refinement criterion
 * like IterativeGridGeneratorRitterNovak, but adapted for the use with
 * fuzzy intervals.
 * Caution: This class uses HashRefinementMultiple, so it generates grids
 * that don't meet the "hierarchical ancestors" requirement!
 *
 * Literature: Julian Valentin. B-Splines for Sparse Grids: Algorithms
 * and Application to Higher-Dimensional Optimization. PhD thesis,
 * University of Stuttgart, IPVS, 2019.
 *
 * @see IterativeGridGeneratorRitterNovak
 * @see HashRefinementMultiple
 */
class IterativeGridGeneratorFuzzyRitterNovak : public IterativeGridGeneratorRitterNovak {
 public:
  static const size_t DEFAULT_NUMBER_OF_ALPHA_SEGMENTS = 10;

  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param f                       objective function
   * @param grid                    grid (should be empty)
   * @param N                       maximal number of grid points
   * @param xFuzzy                  fuzzy input intervals
   * @param numberOfAlphaSegments   number of alpha segments (i.e., number of
   *                                inner Ritter-Novak generations)
   * @param adaptivity              adaptivity between 0 and 1
   * @param initialLevel            level of initial regular sparse grid
   * @param maxLevel                maximal level of grid points
   * @param powMethod               exponentiation method
   *                                (fastPow is faster than std::pow,
   *                                but only approximative)
   */
  IterativeGridGeneratorFuzzyRitterNovak(
      base::ScalarFunction& f, base::Grid& grid, size_t N,
      const std::vector<FuzzyInterval*>& xFuzzy,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS,
      double adaptivity = DEFAULT_ADAPTIVITY,
      base::level_t initialLevel = DEFAULT_INITIAL_LEVEL,
      base::level_t maxLevel = DEFAULT_MAX_LEVEL,
      PowMethod powMethod = PowMethod::STD_POW);

  ~IterativeGridGeneratorFuzzyRitterNovak() override;

  /**
   * Destructor.
   */
  bool generate() override;

  /**
   * @return  fuzzy input intervals
   */
  const std::vector<FuzzyInterval*>& getXFuzzy() const;

  /**
   * @param xFuzzy  fuzzy input intervals
   */
  void setXFuzzy(const std::vector<FuzzyInterval*>& xFuzzy);

  /**
   * @return  number of alpha segments (i.e., number of inner Ritter-Novak
   *          generations)
   */
  size_t getNumberOfAlphaSegments() const;

  /**
   * @param numberOfAlphaSegments   number of alpha segments (i.e., number of
   *                                inner Ritter-Novak generations)
   */
  void setNumberOfAlphaSegments(size_t numberOfAlphaSegments);

 protected:
  /// fuzzy input intervals
  std::vector<FuzzyInterval*> xFuzzy;
  /// number of alpha segments (i.e., number of inner Ritter-Novak generations)
  size_t numberOfAlphaSegments;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORFUZZYRITTERNOVAK_HPP */
