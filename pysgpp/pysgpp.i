// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
%module(directors="1", moduleimport="from . import _pysgpp_swig") pysgpp_swig
// %feature("autodoc", "2");
// %feature("docstring");

%include "base/src/sgpp/globaldef.hpp"

%include "stl.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_pair.i"
%include "std_complex.i"
%include "std_map.i"
%include "std_shared_ptr.i"
%include "carrays.i"
%include "cpointer.i"
%include "typemaps.i"
%include "stdint.i"
%include "exception.i"

%include "carrays.i"
%array_class(unsigned int, unsignedIntArray);
%array_class(bool,BoolArray);
%array_class(int, IntArray);

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
#ifdef _OPENMP
#include <omp.h>
#endif
%}

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}



%{
#ifdef SG_BASE
#include "sgpp/globaldef.hpp"
#include "sgpp_base.hpp"
#endif
#ifdef SG_PDE
#include "sgpp_pde.hpp"
#endif
#ifdef SG_FINANCE
#include "sgpp_finance.hpp"
#endif
#ifdef SG_SOLVER
#include "sgpp_solver.hpp"
#endif
#ifdef SG_PARALLEL
#include "sgpp_parallel.hpp"
#endif
#ifdef SG_QUADRATURE
#include "sgpp_quadrature.hpp"
#endif
#ifdef SG_OPTIMIZATION
#include "sgpp_optimization.hpp"
#endif
#ifdef SG_DATADRIVEN
#include "sgpp_datadriven.hpp"
#endif
#ifdef SG_COMBIGRID
#include "sgpp_combigrid.hpp"
#endif
#ifdef SG_MISC
#include "sgpp_misc.hpp"
#endif
%}


#ifdef SG_BASE
#ifdef PYDOC
%include "base_doc.i"
#endif

%include "base/build/pysgpp/base.i"
#endif

#ifdef SG_PDE
#ifdef PYDOC
%include "pde_doc.i"
#endif

%include "pde/build/pysgpp/pde.i"
#endif

#ifdef SG_FINANCE
#ifdef PYDOC
%include "finance_doc.i"
#endif

%include "finance/build/pysgpp/finance.i"
#endif

#ifdef SG_SOLVER
#ifdef PYDOC
%include "solver_doc.i"
#endif

%include "solver/build/pysgpp/solver.i"
#endif

#ifdef SG_QUADRATURE
#ifdef PYDOC
%include "quadrature_doc.i"
#endif

%include "quadrature/build/pysgpp/quadrature.i"
#endif

#ifdef SG_OPTIMIZATION
#ifdef PYDOC
%include "optimization_doc.i"
#endif

%include "optimization/build/pysgpp/optimization.i"
#endif

#ifdef SG_DATADRIVEN
#ifdef PYDOC
%include "datadriven_doc.i"
#endif

%include "datadriven/build/pysgpp/datadriven.i"
#endif

#ifdef SG_COMBIGRID
#ifdef PYDOC
%include "combigrid_doc.i"
#endif

%include "combigrid/build/pysgpp/combigrid.i"
#endif
