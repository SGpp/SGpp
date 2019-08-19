// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLUnified/OperationMultiEvalStreamingModOCLUnified.hpp>

namespace sgpp {
namespace datadriven {

base::OperationMultipleEval* createStreamingModOCLUnifiedConfigured(
    base::Grid& grid, base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration);

}  // namespace datadriven
}  // namespace sgpp
