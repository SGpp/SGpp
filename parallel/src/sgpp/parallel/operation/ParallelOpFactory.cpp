// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USEMIC
#include <sgpp/parallel/datadriven/basis/common/mic/MICKernel.hpp>
#include <sgpp/parallel/datadriven/basis/common/mic/MICCPUHybridKernel.hpp>
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif
#include <sgpp/parallel/datadriven/basis/linear/noboundary/operation/impl/MICLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/MICModLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/MICModLinearMask.hpp>
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif
#endif

#include <sgpp/parallel/operation/ParallelOpFactory.hpp>

#include <cstring>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalIterative.hpp>
#include <sgpp/parallel/datadriven/basis/linear/noboundary/operation/impl/X86SimdLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/X86SimdModLinearMask.hpp>

#include <sgpp/parallel/datadriven/basis/common/CPUKernel.hpp>
#ifdef USEOCL
#include <sgpp/parallel/datadriven/basis/common/ocl/OCLKernel.hpp>
#include <sgpp/parallel/datadriven/basis/common/ocl/OCLKernelImpl.hpp>
#include <sgpp/parallel/datadriven/basis/common/ocl/OCLCPUHybridKernel.hpp>
#include <sgpp/parallel/datadriven/basis/linear/noboundary/operation/impl/OCLLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/OCLModLinear.hpp>
#include <sgpp/parallel/datadriven/basis/modlinear/operation/impl/OCLModLinearMask.hpp>
#endif

#ifdef USEARBB
#include <sgpp/parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeArBBLinear.hpp>
#endif

#ifdef USE_MPI
#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLaplaceVectorizedLinear.hpp>
#include <sgpp/parallel/pde/basis/linear/boundary/operation/OperationLaplaceVectorizedLinearBoundary.hpp>
#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLTwoDotProductVectorizedLinear.hpp>
#include <sgpp/parallel/pde/basis/linear/boundary/operation/OperationLTwoDotProductVectorizedLinearBoundary.hpp>
//combined operator
#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLTwoDotLaplaceVectorizedLinear.hpp>
#include <sgpp/parallel/pde/basis/linear/boundary/operation/OperationLTwoDotLaplaceVectorizedLinearBoundary.hpp>
#ifdef USEOCL
#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLaplaceVectorizedLinearOCL.hpp>
#include <sgpp/parallel/pde/basis/linear/boundary/operation/OperationLaplaceVectorizedLinearBoundaryOCL.hpp>
#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLTwoDotProductVectorizedLinearOCL.hpp>
#include <sgpp/parallel/pde/basis/linear/boundary/operation/OperationLTwoDotProductVectorizedLinearBoundaryOCL.hpp>
//combined operator
#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLTwoDotLaplaceVectorizedLinearOCL.hpp>
#include <sgpp/parallel/pde/basis/linear/boundary/operation/OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL.hpp>
#endif
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace op_factory {

    parallel::OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(base::Grid& grid, const parallel::VectorizationType vecType, base::DataMatrix* dataset,
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
          return new parallel::OperationMultipleEvalIterative<parallel::MICKernel<parallel::MICLinear> > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterative<parallel::MICCPUHybridKernel<parallel::X86SimdLinear, parallel::MICLinear> > (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#endif
#endif
        else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      }

      else if (strcmp(grid.getType(), "linearL0Boundary") == 0
               || strcmp(grid.getType(), "linearBoundary") == 0) {
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
          return new parallel::OperationMultipleEvalIterative<parallel::MICKernel<parallel::MICLinear> >
                 (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
        }

#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          return new parallel::OperationMultipleEvalIterative<parallel::MICCPUHybridKernel<parallel::X86SimdLinear, parallel::MICLinear> >
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
            return new parallel::OperationMultipleEvalIterative<parallel::CPUKernel<parallel::X86SimdModLinear> >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterative<parallel::CPUKernel<parallel::X86SimdModLinearMask> >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
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
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterative < parallel::MICKernel < parallel::MICModLinear > >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterative < parallel::MICKernel < parallel::MICModLinearMask > >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else {
            throw base::factory_exception("ParallelOpFactory: SGPP_MODLINEAR_EVAL must be 'mask' or 'orig'.");
          }
        }

#ifdef __INTEL_OFFLOAD // Hybrid CPU MIC Mode only makes sense in offload mode
        else if (vecType == parallel::Hybrid_X86SIMD_MIC) {
          if (strcmp(modlinear_mode, "orig") == 0) {
            return new parallel::OperationMultipleEvalIterative < parallel::MICCPUHybridKernel < parallel::X86SimdModLinear, parallel::MICModLinear > >
                   (grid.getStorage(), dataset, gridFrom, gridTo, datasetFrom, datasetTo);
          } else if (strcmp(modlinear_mode, "mask") == 0) {
            return new parallel::OperationMultipleEvalIterative < parallel::MICCPUHybridKernel < parallel::X86SimdModLinearMask, parallel::MICModLinearMask > >
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
        throw base::factory_exception("ParallelOpFactory: OperationMultipleEvalVectorized is not implemented for this grid type.");
      }
    }

#ifdef USE_MPI
    base::OperationMatrix* createOperationLTwoDotProductVectorized(base::Grid& grid, const parallel::VectorizationType& vecType) {
      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLTwoDotProductVectorizedLinear(grid.getStorage());
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLTwoDotProductVectorizedLinearOCL(grid.getStorage());
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLTwoDotProductVectorizedLinearBoundary(grid.getStorage());
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLTwoDotProductVectorizedLinearBoundaryOCL(grid.getStorage());
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else {
        throw base::factory_exception("OperationLTwoDotProductVectorized is not implemented for this grid type.");
      }
    }

    parallel::OperationParabolicPDEMatrixCombined* createOperationLTwoDotLaplaceVectorized(base::Grid& grid, SGPP::base::DataVector& lambda, const parallel::VectorizationType& vecType) {
      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLTwoDotLaplaceVectorizedLinear(grid.getStorage(), lambda);
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLTwoDotLaplaceVectorizedLinearOCL(grid.getStorage(), lambda);
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLTwoDotLaplaceVectorizedLinearBoundary(grid.getStorage(), lambda);
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL(grid.getStorage(), lambda);
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else {
        throw base::factory_exception("OperationLTwoDotLaplaceVectorized is not implemented for this grid type.");
      }
    }

    parallel::OperationParabolicPDEMatrixCombined* createOperationLTwoDotLaplaceVectorized(base::Grid& grid, const parallel::VectorizationType& vecType) {
      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLTwoDotLaplaceVectorizedLinear(grid.getStorage());
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLTwoDotLaplaceVectorizedLinearOCL(grid.getStorage());
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLTwoDotLaplaceVectorizedLinearBoundary(grid.getStorage());
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLTwoDotLaplaceVectorizedLinearBoundaryOCL(grid.getStorage());
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else {
        throw base::factory_exception("OperationLTwoDotLaplaceVectorized is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplaceVectorized(base::Grid& grid, const parallel::VectorizationType& vecType) {
      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLaplaceVectorizedLinear(grid.getStorage());
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLaplaceVectorizedLinearOCL(grid.getStorage());
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLaplaceVectorizedLinearBoundary(grid.getStorage());
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLaplaceVectorizedLinearBoundaryOCL(grid.getStorage());
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else {
        throw base::factory_exception("ParallelOpFactory: OperationLaplaceVectorized is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplaceVectorized(base::Grid& grid, SGPP::base::DataVector& lambda, const parallel::VectorizationType& vecType) {
      if (strcmp(grid.getType(), "linear") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLaplaceVectorizedLinear(grid.getStorage(), lambda);
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLaplaceVectorizedLinearOCL(grid.getStorage(), lambda);
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        if (vecType == parallel::X86SIMD) {
          return new parallel::OperationLaplaceVectorizedLinearBoundary(grid.getStorage(), lambda);
#ifdef USEOCL
        } else if (vecType == parallel::OpenCL) {
          return new parallel::OperationLaplaceVectorizedLinearBoundaryOCL(grid.getStorage(), lambda);
#endif
        } else {
          throw base::factory_exception("Unsupported vectorization type");
        }
      } else {
        throw base::factory_exception("OperationLaplaceVectorized is not implemented for this grid type.");
      }
    }
#endif
  }
}
