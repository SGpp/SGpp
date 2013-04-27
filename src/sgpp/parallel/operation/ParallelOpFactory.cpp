/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
// @author Alexander Heinecke (alexander.heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

//#define USEOCL_NONMASK_MODLINEAR

#include "parallel/operation/ParallelOpFactory.hpp"

#include <cstring>

#include "base/exception/factory_exception.hpp"

#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterative.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSP.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinearMask.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinearMask.hpp"

#ifdef USEARBB
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeArBBLinear.hpp"
#endif
#ifdef USEOCL
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeOCLLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeHybridX86SimdOCLLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeOCLModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeHybridX86SimdOCLModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeOCLModMaskLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeHybridX86SimdOCLModMaskLinear.hpp"
#endif
#ifdef USEMIC
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeMICLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeHybridX86SimdMICLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeMICModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeHybridX86SimdMICModLinear.hpp"
#endif

#ifdef USEARBB
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPArBBLinear.hpp"
#endif
#ifdef USEOCL
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPOCLLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPOCLModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPHybridX86SimdOCLModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPOCLModMaskLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear.hpp"
#endif
#ifdef USEMIC
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPMICLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPHybridX86SimdMICLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPMICModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPHybridX86SimdMICModLinear.hpp"
#endif

#include "parallel/pde/basis/linear/noboundary/operation/OperationLaplaceVectorizedLinear.hpp"
#include "parallel/pde/basis/linear/boundary/operation/OperationLaplaceVectorizedLinearBoundary.hpp"

namespace sg {

  namespace op_factory {


    parallel::OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(base::Grid& grid, const parallel::VectorizationType& vecType, base::DataMatrix* dataset,
        int gridFrom, int gridTo, int datasetFrom, int datasetTo) {
      // handle default upper boundaries
      if (gridTo == -1) gridTo = (int)(grid.getStorage()->size());

      if (datasetTo == -1) datasetTo = (int)(dataset->getNcols());

      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationMultipleEvalIterative<parallel::X86SimdLinear>(grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return new parallel::OperationMultipleEvalIterativeOCLLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return new parallel::OperationMultipleEvalIterativeHybridX86SimdOCLLinear(grid.getStorage(), dataset);
        }

#endif
#ifdef USEARBB
        else if (vecType == parallel::ArBB) {
          return new parallel::OperationMultipleEvalIterativeArBBLinear(grid.getStorage(), dataset);
        }

#endif
#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          return new parallel::OperationMultipleEvalIterativeMICLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterativeHybridX86SimdMICLinear(grid.getStorage(), dataset);
        }

#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      }


      else if (strcmp(grid.getType(), "linearBoundary") == 0
               || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationMultipleEvalIterative<parallel::X86SimdLinear>(
                   grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return new parallel::OperationMultipleEvalIterativeOCLLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return new parallel::OperationMultipleEvalIterativeHybridX86SimdOCLLinear(grid.getStorage(), dataset);
        }

#endif
#ifdef USEARBB
        else if (vecType == parallel::ArBB) {
          return new parallel::OperationMultipleEvalIterativeArBBLinear(grid.getStorage(), dataset);
        }

#endif
#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          return new parallel::OperationMultipleEvalIterativeMICLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterativeHybridX86SimdMICLinear(grid.getStorage(), dataset);
        }

#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else if (strcmp(grid.getType(), "modlinear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationMultipleEvalIterative<parallel::X86SimdModLinear>(grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
#ifdef USEOCL_NONMASK_MODLINEAR
          return new parallel::OperationMultipleEvalIterativeOCLModLinear(grid.getStorage(), dataset);
#else
          return new parallel::OperationMultipleEvalIterativeOCLModMaskLinear(grid.getStorage(), dataset);
#endif
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
#ifdef USEOCL_NONMASK_MODLINEAR
          return new parallel::OperationMultipleEvalIterativeHybridX86SimdOCLModLinear(grid.getStorage(), dataset);
#else
          return new parallel::OperationMultipleEvalIterativeHybridX86SimdOCLModMaskLinear(grid.getStorage(), dataset);
#endif
        }

#endif
#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          return new parallel::OperationMultipleEvalIterativeMICModLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterativeHybridX86SimdMICModLinear(grid.getStorage(), dataset);
        }

#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else {
        throw base::factory_exception("OperationMultipleEvalVectorized is not implemented for this grid type.");
      }
    }

    parallel::OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(base::Grid& grid, const parallel::VectorizationType& vecType, base::DataMatrixSP* dataset,
        int gridFrom, int gridTo, int datasetFrom, int datasetTo) {
      // handle default upper boundaries
      if (gridTo == -1) gridTo = (int)(grid.getStorage()->size());

      if (datasetTo == -1) datasetTo = (int)(dataset->getNcols());

      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPX86SimdLinear>(
                   grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return new parallel::OperationMultipleEvalIterativeSPOCLLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear(grid.getStorage(), dataset);
        }

#endif
#ifdef USEARBB
        else if (vecType == parallel::ArBB) {
          return new parallel::OperationMultipleEvalIterativeSPArBBLinear(grid.getStorage(), dataset);
        }

#endif
#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          return new parallel::OperationMultipleEvalIterativeSPMICLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdMICLinear(grid.getStorage(), dataset);
        }

#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else if (strcmp(grid.getType(), "linearBoundary") == 0
                 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPX86SimdLinear>(
                   grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return new parallel::OperationMultipleEvalIterativeSPOCLLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear(grid.getStorage(), dataset);
        }

#endif
#ifdef USEARBB
        else if (vecType == parallel::ArBB) {
          return new parallel::OperationMultipleEvalIterativeSPArBBLinear(grid.getStorage(), dataset);
        }

#endif
#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          return new parallel::OperationMultipleEvalIterativeSPMICLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdMICLinear(grid.getStorage(), dataset);
        }

#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else if (strcmp(grid.getType(), "modlinear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPX86SimdModLinear>(grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
#ifdef USEOCL_NONMASK_MODLINEAR
          return new parallel::OperationMultipleEvalIterativeSPOCLModLinear(grid.getStorage(), dataset);
#else
          return new parallel::OperationMultipleEvalIterativeSPOCLModMaskLinear(grid.getStorage(), dataset);
#endif
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
#ifdef USEOCL_NONMASK_MODLINEAR
          return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdOCLModLinear(grid.getStorage(), dataset);
#else
          return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdOCLModMaskLinear(grid.getStorage(), dataset);
#endif
        }

#endif
#ifdef USEMIC
        else if (vecType == parallel::MIC) {
          return new parallel::OperationMultipleEvalIterativeSPMICModLinear(grid.getStorage(), dataset);
        } else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdMICModLinear(grid.getStorage(), dataset);
        }

#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else {
        throw base::factory_exception("OperationMultipleEvalVectorizedSP is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplaceVectorized(base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0) {
        return new parallel::OperationLaplaceVectorizedLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new parallel::OperationLaplaceVectorizedLinearBoundary(grid.getStorage());
      } else {
        throw base::factory_exception("OperationLaplaceVectorized is not implemented for this grid type.");
      }
    }
  }
}

