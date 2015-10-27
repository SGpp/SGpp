#pragma once

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class TunableParameter {
private:
    std::string name;
    std::vector<std::string> values;
public:
    TunableParameter(std::string name, std::vector<std::string> values): name(name), values(values) {

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
