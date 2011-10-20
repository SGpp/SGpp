/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "parallel/operation/ParallelOpFactory.hpp"

#include <cstring>

#include "exception/factory_exception.hpp"

#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeX86SimdLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeArBBLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeHybridX86SimdOCLLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalIterativeX86SimdModLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalIterativeOCLModLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalIterativeHybridSSEOCLModLinear.hpp"

#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPX86SimdLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear.hpp"
#include "basis/linear/noboundary/operation/datadriven/OperationMultipleEvalIterativeSPArBBLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalIterativeSPX86SimdModLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalIterativeSPOCLModLinear.hpp"
#include "basis/modlinear/operation/datadriven/OperationMultipleEvalIterativeSPHybridSSEOCLModLinear.hpp"

namespace sg
{

namespace op_factory
{

  
  base::OperationMultipleEvalVectorized* createOperationMultipleEvalVectorized(base::Grid& grid, const std::string& vecType, base::DataMatrix* dataset)
  {

    if(strcmp(grid.getType(), "linear") == 0)
      {
        if (vecType == "X86SIMD")
          {
            return new parallel::OperationMultipleEvalIterativeX86SimdLinear(grid.getStorage(), dataset);
          }
#ifdef USEOCL
        else if (vecType == "OCL")
          {
            return new parallel::OperationMultipleEvalIterativeOCLLinear(grid.getStorage(), dataset);
          }
        else if (vecType == "HYBRID_SSE_OCL")
          {
            return new parallel::OperationMultipleEvalIterativeHybridX86SimdOCLLinear(grid.getStorage(), dataset);
          }
#endif
#ifdef USEARBB
        else if (vecType == "ArBB")
          {
            return new parallel::OperationMultipleEvalIterativeArBBLinear(grid.getStorage(), dataset);
          }
#endif
        else
          {
            throw base::factory_exception("Unsupported vectorization type");
          }
      }


    else if(strcmp(grid.getType(), "linearBoundary") == 0
            || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        if (vecType == "X86SIMD")
          {
            return new parallel::OperationMultipleEvalIterativeX86SimdLinear(grid.getStorage(), dataset);
          }
#ifdef USEOCL
        else if (vecType == "OCL")
          {
            return new parallel::OperationMultipleEvalIterativeOCLLinear(grid.getStorage(), dataset);
          }
        else if (vecType == "HYBRID_SSE_OCL")
          {
            return new parallel::OperationMultipleEvalIterativeHybridX86SimdOCLLinear(grid.getStorage(), dataset);
          }
#endif
#ifdef USEARBB
        else if (vecType == "ArBB")
          {
            return new parallel::OperationMultipleEvalIterativeArBBLinear(grid.getStorage(), dataset);
          }
#endif
        else
          {
            throw base::factory_exception("Unsupported vectorization type");
          }
      }
    else if(strcmp(grid.getType(), "modlinear") == 0)
      {
        if (vecType == "X86SIMD")
          {
            return new parallel::OperationMultipleEvalIterativeX86SimdModLinear(grid.getStorage(), dataset);
          }
#ifdef USEOCL
        else if (vecType == "OCL")
          {
            return new parallel::OperationMultipleEvalIterativeOCLModLinear(grid.getStorage(), dataset);
          }
        else if (vecType == "HYBRID_SSE_OCL")
          {
            return new parallel::OperationMultipleEvalIterativeHybridSSEOCLModLinear(grid.getStorage(), dataset);
          }
#endif
        else
          {
            throw base::factory_exception("Unsupported vectorization type");
          }
      }
    else
      {
        throw base::factory_exception("OperationMultipleEvalVectorized is not implemented for this grid type.");
      }
  }

  base::OperationMultipleEvalVectorizedSP* createOperationMultipleEvalVectorizedSP(base::Grid& grid, const std::string& vecType, base::DataMatrixSP* dataset)
  {
    if(strcmp(grid.getType(), "linear") == 0)
      {
        if (vecType == "X86SIMD")
          {
            return new parallel::OperationMultipleEvalIterativeSPX86SimdLinear(grid.getStorage(), dataset);
          }
#ifdef USEOCL
        else if (vecType == "OCL")
          {
            return new parallel::OperationMultipleEvalIterativeSPOCLLinear(grid.getStorage(), dataset);
          }
        else if (vecType == "HYBRID_SSE_OCL")
          {
            return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear(grid.getStorage(), dataset);
          }
#endif
#ifdef USEARBB
        else if (vecType == "ArBB")
          {
            return new parallel::OperationMultipleEvalIterativeSPArBBLinear(grid.getStorage(), dataset);
          }
#endif
        else
          {
            throw base::factory_exception("Unsupported vectorization type");
          }
      }
    else if(strcmp(grid.getType(), "linearBoundary") == 0
            || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0)
      {
        if (vecType == "X86SIMD")
          {
            return new parallel::OperationMultipleEvalIterativeSPX86SimdLinear(grid.getStorage(), dataset);
          }
#ifdef USEOCL
        else if (vecType == "OCL")
          {
            return new parallel::OperationMultipleEvalIterativeSPOCLLinear(grid.getStorage(), dataset);
          }
        else if (vecType == "HYBRID_SSE_OCL")
          {
            return new parallel::OperationMultipleEvalIterativeSPHybridX86SimdOCLLinear(grid.getStorage(), dataset);
          }
#endif
#ifdef USEARBB
        else if (vecType == "ArBB")
          {
            return new parallel::OperationMultipleEvalIterativeSPArBBLinear(grid.getStorage(), dataset);
          }
#endif
        else
          {
            throw base::factory_exception("Unsupported vectorization type");
          }
      }
    else if(strcmp(grid.getType(), "modlinear") == 0)
      {
        if (vecType == "X86SIMD")
          {
            return new parallel::OperationMultipleEvalIterativeSPX86SimdModLinear(grid.getStorage(), dataset);
          }
#ifdef USEOCL
        else if (vecType == "OCL")
          {
            return new parallel::OperationMultipleEvalIterativeSPOCLModLinear(grid.getStorage(), dataset);
          }
        else if (vecType == "HYBRID_SSE_OCL")
          {
            return new parallel::OperationMultipleEvalIterativeSPHybridSSEOCLModLinear(grid.getStorage(), dataset);
          }
#endif
        else
          {
            throw base::factory_exception("Unsupported vectorization type");
          }
      }
    else
      {
        throw base::factory_exception("OperationMultipleEvalVectorizedSP is not implemented for this grid type.");
      }
  }

}
}

