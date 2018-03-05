// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

%module(directors="1") matsgpp

%include "base/src/sgpp/globaldef.hpp"

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_string.i"
%include "stdint.i"

%include "typemaps.i"

%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
#include <omp.h>
%}

%{
#include "sgpp_base.hpp"
#ifdef SG_PDE
#include "sgpp_pde.hpp"
#endif
#ifdef SG_FINANCE
#include "sgpp_finance.hpp"
#endif
#ifdef SG_SOLVER
#include "sgpp_solver.hpp"
#endif
#ifdef SG_DATADRIVEN
#include "sgpp_datadriven.hpp"
#endif
#ifdef SG_PARALLEL
#include "sgpp_parallel.hpp"
#endif
#ifdef SG_COMBIGRID
#include "sgpp_combigrid.hpp"
#endif
#ifdef SG_QUADRATURE
#include "sgpp_quadrature.hpp"
#endif
#ifdef SG_OPTIMIZATION
#include "sgpp_optimization.hpp"
#endif
#ifdef SG_MISC
#include "sgpp_misc.hpp"
#endif
%}

%include "base/build/matsgpp/base.i"

#ifdef SG_PDE
%include "pde/build/matsgpp/pde.i"
#endif

#ifdef SG_FINANCE
%include "finance/build/matsgpp/finance.i"
#endif

#ifdef SG_SOLVER
%include "solver/build/matsgpp/solver.i"
#endif

#ifdef SG_DATADRIVEN
%include "datadriven/build/matsgpp/datadriven.i"
#endif

#ifdef SG_QUADRATURE
%include "quadrature/build/matsgpp/quadrature.i"
#endif

#ifdef SG_COMBIGRID
%include "combigrid/build/matsgpp/combigrid.i"
#endif

#ifdef SG_OPTIMIZATION
%include "optimization/build/matsgpp/optimization.i"
#endif
