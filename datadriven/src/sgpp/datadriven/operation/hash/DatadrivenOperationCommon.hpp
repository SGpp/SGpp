// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/OperationConfiguration.hpp>
#include <sgpp/globaldef.hpp>

#include <memory>
#include <string>

namespace sgpp {
namespace datadriven {

enum class OperationMultipleEvalType {
  DEFAULT,
  STREAMING,
  SUBSPACELINEAR,
  ADAPTIVE,
  MORTONORDER,
  SCALAPACK
};

enum class OperationMultipleEvalSubType {
  DEFAULT,
  SIMPLE,
  COMBINED,
  OCL,
  OCLFASTMP,
  OCLMP,
  OCLMASKMP,
  OCLOPT,
  OCLUNIFIED,
  CUDA
};

enum class OperationMultipleEvalMPIType { NONE, MASTERSLAVE, HPX };

class OperationMultipleEvalConfiguration {
 private:
  OperationMultipleEvalType type = OperationMultipleEvalType::DEFAULT;
  OperationMultipleEvalSubType subType = OperationMultipleEvalSubType::DEFAULT;
  OperationMultipleEvalMPIType mpiType = OperationMultipleEvalMPIType::NONE;

  std::shared_ptr<base::OperationConfiguration> parameters;

  // optional - can be set for easier reporting
  std::string name;

 public:
  OperationMultipleEvalConfiguration(
      OperationMultipleEvalType type = OperationMultipleEvalType::DEFAULT,
      OperationMultipleEvalSubType subType = OperationMultipleEvalSubType::DEFAULT,
      OperationMultipleEvalMPIType mpiType = OperationMultipleEvalMPIType::NONE,
      std::string name = "unnamed") {
    this->type = type;
    this->subType = subType;
    this->mpiType = mpiType;
    this->name = name;
  }

  OperationMultipleEvalConfiguration(OperationMultipleEvalType type,
                                     OperationMultipleEvalSubType subType,
                                     base::OperationConfiguration& parameters,
                                     std::string name = "unnamed") {
    this->type = type;
    this->subType = subType;
    this->name = name;
    this->parameters = std::shared_ptr<base::OperationConfiguration>(parameters.clone());
  }

  OperationMultipleEvalConfiguration(OperationMultipleEvalType type,
                                     OperationMultipleEvalSubType subType,
                                     OperationMultipleEvalMPIType mpiType,
                                     base::OperationConfiguration& parameters,
                                     std::string name = "unnamed") {
    this->type = type;
    this->subType = subType;
    this->mpiType = mpiType;
    this->name = name;
    this->parameters = std::shared_ptr<base::OperationConfiguration>(parameters.clone());
  }

  OperationMultipleEvalMPIType getMPIType() { return this->mpiType; }

  OperationMultipleEvalType getType() { return this->type; }

  OperationMultipleEvalSubType getSubType() { return this->subType; }

  std::shared_ptr<base::OperationConfiguration> getParameters() { return this->parameters; }

  const std::string& getName() const { return this->name; }
};
}  // namespace datadriven
}  // namespace sgpp
