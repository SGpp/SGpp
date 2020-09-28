// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TYPESSOLVER_HPP
#define TYPESSOLVER_HPP

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace solver {

/**
 * enum to address different SLE solvers in a standardized way
 */
enum class SLESolverType { CG, BiCGSTAB, FISTA };

struct SLESolverConfiguration {
  sgpp::solver::SLESolverType type_ = SLESolverType::CG;
  double eps_ = 1e-12;
  size_t maxIterations_ = 100;
  double threshold_ = 1e-12;
  bool verbose_ = false;
};

struct SLESolverSPConfiguration {
  sgpp::solver::SLESolverType type_;
  float eps_;
  size_t maxIterations_;
  float threshold_;
};
}  // namespace solver
}  // namespace sgpp

#endif /* TYPESSOLVER_HPP */
