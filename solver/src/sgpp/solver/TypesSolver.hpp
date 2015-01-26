/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TYPESSOLVER_HPP
#define TYPESSOLVER_HPP

#include <cstddef>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace solver {

    /**
     * enum to address different SLE solvers in a standardized way
     */
    enum SLESolverType {
      CG,
      BiCGSTAB
    };

    struct SLESolverConfiguration {
      SGPP::solver::SLESolverType type_;
      double eps_;
      size_t maxIterations_;
      double threshold_;
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
