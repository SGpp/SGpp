// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingOCLMultiPlatformConfigured(
    base::Grid& grid, base::DataMatrix& dataset,
    SGPP::datadriven::OperationMultipleEvalConfiguration& configuration);

}  // namespace datadriven
}  // namespace SGPP
