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

#include "parallel/datadriven/operation/OperationMultipleEvalIterative.hpp"
#include "parallel/datadriven/operation/OperationMultipleEvalIterativeSP.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinear.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinearMask.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinearMask.hpp"

#include "parallel/datadriven/basis/common/CPUKernel.hpp"
#include "parallel/datadriven/basis/common/SPCPUKernel.hpp"
#ifdef USEOCL
#include "parallel/datadriven/basis/common/OCLKernel.hpp"
#include "parallel/datadriven/basis/common/OCLKernelImpl.hpp"
#include "parallel/datadriven/basis/common/OCLCPUHybridKernel.hpp"
#include "parallel/datadriven/basis/common/SPOCLKernel.hpp"
#include "parallel/datadriven/basis/common/SPOCLKernelImpl.hpp"
#include "parallel/datadriven/basis/common/SPOCLCPUHybridKernel.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/OCLLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/OCLModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/OCLModLinearMask.hpp"
#endif
#ifdef USEARBB
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeArBBLinear.hpp"
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
        size_t gridFrom, size_t gridTo, size_t datasetFrom, size_t datasetTo) {
      // handle default upper boundaries
      if (gridTo == std::numeric_limits<size_t>::max()) {
        gridTo = grid.getStorage()->size();
      }

      if (datasetTo == std::numeric_limits<size_t>::max()) {
        datasetTo = dataset->getNcols();
      }

      // get env var
      const char* modlinear_mode = getenv("SGPP_MODLINEAR_EVAL");

      if (modlinear_mode == NULL) {
        modlinear_mode = "mask";
      }

      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationMultipleEvalIterative<parallel::CPUKernel<parallel::X86SimdLinear> >(grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return new parallel::OperationMultipleEvalIterative <
                 parallel::OCLKernel<parallel::OCLLinear<double> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return new parallel::OperationMultipleEvalIterative < parallel::OCLCPUHybridKernel < parallel::X86SimdLinear,
                 parallel::OCLLinear<double> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
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
          return new parallel::OperationMultipleEvalIterative<parallel::CPUKernel<parallel::X86SimdLinear> >(
                   grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return new parallel::OperationMultipleEvalIterative <
                 parallel::OCLKernel < parallel::OCLLinear<double> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return new parallel::OperationMultipleEvalIterative <
                 parallel::OCLCPUHybridKernel < parallel::X86SimdLinear,
                 parallel::OCLLinear<double> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
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
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterative<parallel::CPUKernel<parallel::X86SimdModLinear> >(grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterative<parallel::CPUKernel<parallel::X86SimdModLinearMask> >(grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterative <
                   parallel::OCLKernel < parallel::OCLModLinear<double> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterative <
                   parallel::OCLKernel < parallel::OCLModLinearMask<double> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }

        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterative <
                   parallel::OCLCPUHybridKernel < parallel::X86SimdModLinear,
                   parallel::OCLModLinear<double> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterative <
                   parallel::OCLCPUHybridKernel < parallel::X86SimdModLinearMask,
                   parallel::OCLModLinearMask<double> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
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
        throw base::factory_exception("ParallelOpFactory: OperationMultipleEvalVectorized is not implemented for this grid type.");
      }
    }

    parallel::OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(base::Grid& grid, const parallel::VectorizationType& vecType, base::DataMatrixSP* dataset,
        size_t gridFrom, size_t gridTo, size_t datasetFrom, size_t datasetTo) {
      // handle default upper boundaries
      if (gridTo == std::numeric_limits<size_t>::max()) {
        gridTo = grid.getStorage()->size();
      }

      if (datasetTo == std::numeric_limits<size_t>::max()) {
        datasetTo = dataset->getNcols();
      }

      // get env var
      const char* modlinear_mode = getenv("SGPP_MODLINEAR_EVAL");

      if (modlinear_mode == NULL) {
        modlinear_mode = "mask";
      }

      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPCPUKernel<parallel::SPX86SimdLinear> >(
                   grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return new parallel::OperationMultipleEvalIterativeSP <
                 parallel::SPOCLKernel < parallel::OCLLinear<float> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return new parallel::OperationMultipleEvalIterativeSP <
                 parallel::SPOCLCPUHybridKernel < parallel::SPX86SimdLinear, parallel::OCLLinear<float> >  > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
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
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPCPUKernel<parallel::SPX86SimdLinear> >(
                   grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          return new parallel::OperationMultipleEvalIterativeSP <
                 parallel::SPOCLKernel < parallel::OCLLinear<float>  > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          return new parallel::OperationMultipleEvalIterativeSP <
                 parallel::SPOCLCPUHybridKernel < parallel::SPX86SimdLinear, parallel::OCLLinear<float>  > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
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
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP<parallel::SPCPUKernel<parallel::SPX86SimdModLinear> >(grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP<parallel::SPCPUKernel<parallel::SPX86SimdModLinearMask> >(grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#ifdef USEOCL
        else if (vecType == parallel::OpenCL) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP <
                   parallel::SPOCLKernel < parallel::OCLModLinear<float > > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP <
                   parallel::SPOCLKernel < parallel::OCLModLinearMask<float> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        } else if (vecType == parallel::Hybrid_X86SIMD_OpenCL) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP <
                   parallel::SPOCLCPUHybridKernel < parallel::SPX86SimdModLinear,
                   parallel::OCLModLinear<float> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP <
                   parallel::SPOCLCPUHybridKernel < parallel::SPX86SimdModLinearMask,
                   parallel::OCLModLinearMask<float> > > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
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
        throw base::factory_exception("ParallelOpFactory: OperationMultipleEvalVectorizedSP is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplaceVectorized(base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0) {
        return new parallel::OperationLaplaceVectorizedLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new parallel::OperationLaplaceVectorizedLinearBoundary(grid.getStorage());
      } else {
        throw base::factory_exception("ParallelOpFactory: OperationLaplaceVectorized is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplaceVectorized(base::Grid& grid, sg::base::DataVector& lambda) {
      if (strcmp(grid.getType(), "linear") == 0) {
        return new parallel::OperationLaplaceVectorizedLinear(grid.getStorage(), lambda);
      } else if (strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new parallel::OperationLaplaceVectorizedLinearBoundary(grid.getStorage(), lambda);
      } else {
        throw base::factory_exception("OperationLaplaceVectorized is not implemented for this grid type.");
      }
    }
  }
}

