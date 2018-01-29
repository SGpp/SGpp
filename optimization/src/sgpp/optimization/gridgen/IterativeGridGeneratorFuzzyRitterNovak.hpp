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

class IterativeGridGeneratorFuzzyRitterNovak : public IterativeGridGeneratorRitterNovak {
 public:
  static const size_t DEFAULT_NUMBER_OF_ALPHA_SEGMENTS = 10;

  IterativeGridGeneratorFuzzyRitterNovak(
      ScalarFunction& f, base::Grid& grid, size_t N,
      const std::vector<const FuzzyInterval*>& xFuzzy,
      size_t numberOfAlphaSegments = DEFAULT_NUMBER_OF_ALPHA_SEGMENTS,
      double adaptivity = DEFAULT_ADAPTIVITY,
      base::level_t initialLevel = DEFAULT_INITIAL_LEVEL,
      base::level_t maxLevel = DEFAULT_MAX_LEVEL,
      PowMethod powMethod = STD_POW);

  ~IterativeGridGeneratorFuzzyRitterNovak() override;

  bool generate() override;

  const std::vector<const FuzzyInterval*>& getXFuzzy() const;
  void setXFuzzy(const std::vector<const FuzzyInterval*>& xFuzzy);

  size_t getNumberOfAlphaSegments() const;
  void setNumberOfAlphaSegments(size_t numberOfAlphaSegments);

 protected:
  std::vector<const FuzzyInterval*> xFuzzy;
  size_t numberOfAlphaSegments;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_GRIDGEN_ITERATIVEGRIDGENERATORFUZZYRITTERNOVAK_HPP */
