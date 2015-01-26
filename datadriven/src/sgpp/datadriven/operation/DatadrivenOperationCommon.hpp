#pragma once

#include <string>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

enum class OperationMultipleEvalType {DEFAULT, STREAMING, SUBSPACELINEAR};

enum class OperationMultipleEvalSubType {DEFAULT, SIMPLE, COMBINED};

class OperationMultipleEvalConfiguration {
public:
	OperationMultipleEvalType type = OperationMultipleEvalType::DEFAULT;
	OperationMultipleEvalSubType subType = OperationMultipleEvalSubType::DEFAULT;

	//operational - can be set for easier reporting
	std::string name = "DEFAULT";
};

}
}
