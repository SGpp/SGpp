#pragma once

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

enum class ParameterType {
    ID, UINT, INT, DOUBLE, TEXT, BOOL
};

class TunableParameter {
private:
    std::string name;
    std::vector<std::string> values;
    ParameterType type;
public:
    TunableParameter(std::string name, std::vector<std::string> values, ParameterType type) :
            name(name), values(values), type(type) {

    }

    std::string &getName() {
        return this->name;
    }

    std::vector<std::string> &getValues() {
        return this->values;
    }

};

}
}
