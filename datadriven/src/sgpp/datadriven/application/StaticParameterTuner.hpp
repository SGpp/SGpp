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

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/TunableParameter.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>

namespace SGPP {
namespace datadriven {

class StaticParameterTuner {
private:
    bool collectStatistics;
    bool verbose;
    std::vector<std::pair<SGPP::base::OCLOperationConfiguration, double>> statistics;SGPP::base::OCLOperationConfiguration fixedParameters;
    std::vector<TunableParameter> tunableParameters;

    double evaluateSetup(SGPP::datadriven::LearnerScenario &scenario,
            SGPP::base::OCLOperationConfiguration &currentParameters);
public:
    StaticParameterTuner(bool collectStatistics = false, bool verbose = false);

    //write a file with detailed stats for the optimization
    StaticParameterTuner(std::string tunerFileName, bool collectStatistics = false, bool verbose = false);

    void addFixedParameter(const std::string &name, const std::string &value, const ParameterType type);

    void addParameter(const std::string &name, const std::vector<std::string> &valueRange, const ParameterType type);

    SGPP::base::OCLOperationConfiguration tuneParameters(SGPP::datadriven::LearnerScenario &scenario);

    void writeToFile(const std::string &fileName);

    void readFromFile(const std::string &fileName);

    void writeStatisticsToFile(const std::string &statisticsFileName);

};

}
}

#endif
