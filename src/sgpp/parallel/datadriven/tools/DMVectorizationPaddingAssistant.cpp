/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"

namespace sg
{

namespace parallel
{

size_t DMVectorizationPaddingAssistant::padDataset(sg::base::DataMatrix& dataset, VectorizationType& vecType)
{
    size_t vecWidth = 0;

	if (vecType == X86SIMD)
    {
        vecWidth = 24;
    }
    else if (vecType == OpenCL)
    {
        vecWidth = 128;
    }
    else if (vecType == Hybrid_X86SIMD_OpenCL)
    {
        vecWidth = 128;
    }
    else if (vecType == ArBB)
    {
        vecWidth = 16;
    }
    else if (vecType == MIC)
    {
        vecWidth = 96;
    }
    else if (vecType == Hybrid_X86SIMD_MIC)
    {
        vecWidth = 96;
    }
    else
    {
        throw new sg::base::operation_exception("DMVectorizationPaddingAssistant::padDataset : un-supported vector extension!");
    }

	// Assure that data has a even number of instances -> padding might be needed
    size_t remainder = dataset.getNrows() % vecWidth;
    size_t loopCount = vecWidth - remainder;

    if (loopCount != vecWidth)
    {
        sg::base::DataVector lastRow(dataset.getNcols());
        size_t oldSize = dataset.getNrows();
        dataset.getRow(dataset.getNrows()-1, lastRow);
        dataset.resize(dataset.getNrows()+loopCount);
        for (size_t i = 0; i < loopCount; i++)
        {
            dataset.setRow(oldSize+i, lastRow);
        }
    }

    return dataset.getNrows();
}

size_t DMVectorizationPaddingAssistant::padDataset(sg::base::DataMatrixSP& dataset, VectorizationType vecType)
{
    size_t vecWidth = 0;

	if (vecType == X86SIMD)
    {
        vecWidth = 48;
    }
    else if (vecType == OpenCL)
    {
        vecWidth = 128;
    }
    else if (vecType == Hybrid_X86SIMD_OpenCL)
    {
        vecWidth = 128;
    }
    else if (vecType == ArBB)
    {
        vecWidth = 16;
    }
    else if (vecType == MIC)
    {
        vecWidth = 192;
    }
    else if (vecType == Hybrid_X86SIMD_MIC)
    {
        vecWidth = 192;
    }
    else
    {
        throw new sg::base::operation_exception("DMVectorizationPaddingAssistant::padDataset : un-supported vector extension!");
    }

	// Assure that data has a even number of instances -> padding might be needed
    size_t remainder = dataset.getNrows() % vecWidth;
    size_t loopCount = vecWidth - remainder;

    if (loopCount != vecWidth)
    {
        sg::base::DataVectorSP lastRow(dataset.getNcols());
        size_t oldSize = dataset.getNrows();
        dataset.getRow(dataset.getNrows()-1, lastRow);
        dataset.resize(dataset.getNrows()+loopCount);
        for (size_t i = 0; i < loopCount; i++)
        {
            dataset.setRow(oldSize+i, lastRow);
        }
    }

    return dataset.getNrows();
}

}

}
