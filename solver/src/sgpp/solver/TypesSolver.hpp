// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef TYPESSOLVER_HPP
#define TYPESSOLVER_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace solver {

    /**
     * enum to address different SLE solvers in a standardized way
     */
    enum class SLESolverType {
      CG,
      BiCGSTAB
    };

    struct SLESolverConfiguration {
      SGPP::solver::SLESolverType type_;
      float_t eps_;
      size_t maxIterations_;
      float_t threshold_;
    };

    struct SLESolverSPConfiguration {
      SGPP::solver::SLESolverType type_;
      float eps_;
      size_t maxIterations_;
      float threshold_;
    };

  }

}

#endif /* TYPESSOLVER_HPP */
