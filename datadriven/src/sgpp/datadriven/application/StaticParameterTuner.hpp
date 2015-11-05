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

#if USE_OCL == 1

#include <sgpp/globaldef.hpp>

#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/datadriven/application/TunableParameter.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>

namespace SGPP {
namespace datadriven {
class StaticParameterTuner {
private:
    bool collectStatistics;
    bool verbose;
    std::vector<std::pair<SGPP::base::OCLConfigurationParameters, double>> statistics;
    SGPP::base::OCLConfigurationParameters fixedParameters;
    std::vector<TunableParameter> tunableParameters;

    double evaluateSetup(SGPP::datadriven::LearnerScenario scenario,
    SGPP::base::OCLConfigurationParameters currentParameters);
public:
    StaticParameterTuner(bool collectStatistics = false, bool verbose = false);

    //write a file with detailed stats for the optimization
    StaticParameterTuner(std::string tunerFileName, bool collectStatistics = false, bool verbose = false);

    void addFixedParameter(std::string name, std::string value);

    void addParameter(std::string name, std::vector<std::string> valueRange);

    SGPP::base::OCLConfigurationParameters tuneParameters(SGPP::datadriven::LearnerScenario scenario);

    void writeToFile(std::string fileName);

    void readFromFile(std::string fileName);

    void writeStatisticsToFile(std::string statisticsFileName);

};

}
}

#endif
