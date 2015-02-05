%module(directors="1") pysgpp
%feature("docstring");

%include "../base/src/sgpp/globaldef.hpp"

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_complex.i"
%include "std_map.i"
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
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}



%{

#include "sgpp_base.hpp"
using namespace SGPP::base;
#ifdef SG_DATADRIVEN
#include "sgpp_datadriven.hpp"
using namespace SGPP::datadriven;
#endif
#ifdef SG_PDE
#include "sgpp_pde.hpp"
using namespace SGPP::pde;
#endif
#ifdef SG_FINANCE
#include "sgpp_finance.hpp"
using namespace SGPP::finance;
#endif
#ifdef SG_SOLVER
#include "sgpp_solver.hpp"
using namespace SGPP::solver;
#endif
#ifdef SG_PARALLEL
#include "sgpp_parallel.hpp"
using namespace SGPP::parallel;
#endif
#ifdef SG_COMBIGRID
#include "combigrid.hpp"
using namespace combigrid;
#endif
#ifdef SG_QUADRATURE
#include "sgpp_quadrature.hpp"
using namespace SGPP::quadrature;
#endif
#ifdef SG_OPTIMIZATION
#include "sgpp_optimization.hpp"
#endif
#ifdef SG_MISC
#include "sgpp_misc.hpp"
#endif

using namespace SGPP;

%}


%include "base/build/pysgpp/base.i"

#ifdef SG_DATADRIVEN
%include "datadriven/build/pysgpp/datadriven.i"
#endif

#ifdef SG_PDE
%include "pde/build/pysgpp/pde.i"
#endif

#ifdef SG_FINANCE
%include "finance/build/pysgpp/finance.i"
#endif

#ifdef SG_SOLVER
%include "solver/build/pysgpp/solver.i"
#endif

#ifdef SG_QUADRATURE
%include "quadrature/build/pysgpp/quadrature.i"
#endif

#ifdef SG_COMBIGRID
%include "combigrid/build/pysgpp/combigrid.i"
#endif

#ifdef SG_OPTIMIZATION
%include "optimization/build/pysgpp/optimization.i"
#endif
