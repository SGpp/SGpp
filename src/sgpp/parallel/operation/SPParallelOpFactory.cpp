/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
// @author Alexander Heinecke (alexander.heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifdef USEMIC
#include "parallel/datadriven/basis/common/mic/SPMICKernel.hpp"
#include "parallel/datadriven/basis/common/mic/SPMICCPUHybridKernel.hpp"
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPMICLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPMICModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPMICModLinearMask.hpp"
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
#endif

#include "parallel/operation/ParallelOpFactory.hpp"

#include <cstring>

#include "base/exception/factory_exception.hpp"

#include "parallel/datadriven/operation/OperationMultipleEvalIterativeSP.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/SPX86SimdLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinearMask.hpp"

#include "parallel/datadriven/basis/common/SPCPUKernel.hpp"
#ifdef USEOCL
#include "parallel/datadriven/basis/common/ocl/SPOCLKernel.hpp"
#include "parallel/datadriven/basis/common/ocl/SPOCLKernelImpl.hpp"
#include "parallel/datadriven/basis/common/ocl/SPOCLCPUHybridKernel.hpp"
#include "parallel/datadriven/basis/linear/noboundary/operation/impl/OCLLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/OCLModLinear.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/OCLModLinearMask.hpp"
#endif

#ifdef USEARBB
#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPArBBLinear.hpp"
#endif

namespace sg {

  namespace op_factory {

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
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPMICKernel<parallel::SPMICLinear> > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }
#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPMICCPUHybridKernel<parallel::SPX86SimdLinear, parallel::SPMICLinear> > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }
#endif
#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      }

      else if (strcmp(grid.getType(), "linearBoundary") == 0
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
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPMICKernel<parallel::SPMICLinear> >
                 (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterativeSP<parallel::SPMICCPUHybridKernel<parallel::SPX86SimdLinear, parallel::SPMICLinear> >
                 (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#endif
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
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP < parallel::SPMICKernel < parallel::SPMICModLinear > >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP < parallel::SPMICKernel < parallel::SPMICModLinearMask > >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP < parallel::SPMICCPUHybridKernel < parallel::SPX86SimdModLinear, parallel::SPMICModLinear > >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterativeSP < parallel::SPMICCPUHybridKernel < parallel::SPX86SimdModLinearMask, parallel::SPMICModLinearMask > >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#endif
#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else {
        throw base::factory_exception("ParallelOpFactory: OperationMultipleEvalVectorizedSP is not implemented for this grid type.");
      }
    }
  }
}

