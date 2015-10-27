/*
 * StaticParameterTuner.hpp
 *
 *  Created on: Oct 27, 2015
 *      Author: pfandedd
 */

/*
 * StaticParameterTuner.hpp
 *
 *  Created on: Oct 27, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/datadriven/application/TunableParameter.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>

namespace SGPP {
namespace datadriven {
class StaticParameterTuner {
private:
    bool writeStatistics;
    std::string outFileName;
    SGPP::base::OCLConfigurationParameters fixedParameters;
    std::vector<TunableParameter> tunableParameters;

    double evaluateSetup(SGPP::datadriven::LearnerScenario scenario,
    SGPP::base::OCLConfigurationParameters currentParameters);
public:
    StaticParameterTuner();

    //write a file with detailed stats for the optimization
    StaticParameterTuner(std::string outFileName);

    void addFixedParameter(std::string name, std::string value);

    void addParameter(std::string name, std::vector<std::string> valueRange);

    SGPP::base::OCLConfigurationParameters tuneParameters(SGPP::datadriven::LearnerScenario scenario);

};

}
}

