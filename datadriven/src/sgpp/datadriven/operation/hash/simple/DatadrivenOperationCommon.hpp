// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <string>
#include <memory>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/OperationConfiguration.hpp>

namespace SGPP {
namespace datadriven {

enum class OperationMultipleEvalType {
    DEFAULT, STREAMING, SUBSPACELINEAR, ADAPTIVE
};

enum class OperationMultipleEvalSubType {
    DEFAULT, SIMPLE, COMBINED, OCL, OCLFAST, OCLFASTMULTIPLATFORM, OCLMP, OCLMASK, OCLMASKMP
};

//TODO: remove pointer and create constructor with parameters
class OperationMultipleEvalConfiguration {
private:
    OperationMultipleEvalType type = OperationMultipleEvalType::DEFAULT;
    OperationMultipleEvalSubType subType = OperationMultipleEvalSubType::DEFAULT;
    std::shared_ptr<base::OperationConfiguration> parameters;

    //optional - can be set for easier reporting
    std::string name;
public:

    OperationMultipleEvalConfiguration(OperationMultipleEvalType type = OperationMultipleEvalType::DEFAULT,
            OperationMultipleEvalSubType subType = OperationMultipleEvalSubType::DEFAULT, std::string name = "unnamed") {
        this->type = type;
        this->subType = subType;
        this->name = name;
    }

    OperationMultipleEvalConfiguration(OperationMultipleEvalType type, OperationMultipleEvalSubType subType,
            base::OperationConfiguration &parameters, std::string name = "unnamed") {
        this->type = type;
        this->subType = subType;
        this->name = name;
        this->parameters = std::shared_ptr<base::OperationConfiguration>(parameters.clone());
    }

    OperationMultipleEvalType getType() {
        return this->type;
    }

    OperationMultipleEvalSubType getSubType() {
        return this->subType;
    }

    //TODO: change this to return a reference
    std::shared_ptr<base::OperationConfiguration> getParameters() {
        return this->parameters;
    }

    std::string &getName() {
        return this->name;
    }

};

}
}
