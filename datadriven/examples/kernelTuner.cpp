/*
 * kernelTuner.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: pfandedd
 */
#include <iostream>

#if USE_OCL == 1
#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/StaticParameterTuner.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>

int main(int argc, char **argv) {
/*
    std::string scenarioFileName;
    std::string parameterConfigurationFileName; // <- includes kernel type, subtype
    std::vector<std::string> devices;

    boost::program_options::options_description description("Allowed options");
    description.add_options()("help", "display help")("scenario",
            boost::program_options::value<std::string>(&scenarioFileName),
            "the scenario file to be used (serialized LearnerScenario)")("parameterConfiguration",
            boost::program_options::value<std::string>(&parameterConfigurationFileName),
            "the parameter configuration file to be used (serialized StaticParameterTuner)")("devices",
            boost::program_options::value<std::vector<std::string> >(&devices)->multitoken(),
            "specify comma-separated list of devices or \"all\" for all devices found");

    boost::program_options::variables_map variables_map;

    boost::program_options::parsed_options options = parse_command_line(argc, argv, description);
    boost::program_options::store(options, variables_map);
    boost::program_options::notify(variables_map);

    if (variables_map.count("help")) {
        std::cout << description << std::endl;
        return 0;
    }

    //check whether all files exist
    std::ifstream scenarioFile(scenarioFileName);
    if (!scenarioFile.good()) {
        std::cout << "scenario file not found" << std::endl;
        return 1;
    }
    scenarioFile.close();
    std::ifstream parameterConfigurationFile(parameterConfigurationFileName);
    if (!parameterConfigurationFile.good()) {
        std::cout << "scenario file not found" << std::endl;
        return 1;
    }
    parameterConfigurationFile.close();

    SGPP::datadriven::LearnerScenario scenario(scenarioFileName);

    SGPP::datadriven::StaticParameterTuner staticParameterTuner(parameterConfigurationFileName, true, true);

    for (std::string device : devices) {
        SGPP::base::OCLConfigurationParameters bestParameters = staticParameterTuner.tuneParameters(device, scenario);
        std::string uniqueName = staticParameterTuner.getUniqueName();
        staticParameterTuner.writeStatisticsToFile(uniqueName + ".csv");
        bestParameters.writeToFile("optimal_" + uniqueName + ".cfg");
    }
*/
    return 0;
}
#else
int main(int argc, char **argv) {
    std::cout << "no OpenCL support" << std::endl;
}
#endif

