// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#pragma once

#include <map>

//#include <sgpp/datadriven/operation/OperationMultipleEvalSubspace/simple/X86Simple.hpp>
//#include "X86Baseline.hpp"
//#include "X86Vectorized.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

// list of kernels
enum ComputeKernelType {
	SIMPLE, COMBINED
};

static std::map<ComputeKernelType, std::string> ComputeKernelTypeMap = { { SIMPLE, "SIMPLE" }, { COMBINED, "COMBINED" } };

static std::map<std::string, ComputeKernelType> NameComputeKernelTypeMap = { { "SIMPLE", SIMPLE }, { "COMBINED", COMBINED } };

}
}

//TODO: add vector selection here, basis for auto tuning => use only easily parsable structures
