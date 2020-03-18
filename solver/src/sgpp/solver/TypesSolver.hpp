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
  sgpp::solver::SLESolverType type_;
  double eps_;
  size_t maxIterations_;
  double threshold_;
  bool verbose_;

  /*
  // Debug method to neatly print internal data
  void dumpToStream(std::ostream& stream_out = std::cout) const {
    stream_out << "type: \t\t" << SLESolverTypeParser::toString(type_)
               << std::endl;
    stream_out << "eps: \t\t" << eps_ << std::endl;
    stream_out << "maxIterations: \t" << maxIterations_ << std::endl;
    stream_out << "threshold: \t" << threshold_ << std::endl;
    stream_out << "verbose: \t" << std::boolalpha << verbose_ << std::endl;
  }
  */
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
