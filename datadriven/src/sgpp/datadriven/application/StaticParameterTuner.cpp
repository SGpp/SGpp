#if USE_OCL == 1

#include <limits>
#include <algorithm>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/StaticParameterTuner.hpp>
#include "sgpp/datadriven/application/MetaLearner.hpp"

namespace SGPP {
namespace datadriven {

StaticParameterTuner::StaticParameterTuner(SGPP::base::OCLOperationConfiguration &fixedParameters,
        bool collectStatistics, bool verbose) :
        collectStatistics(collectStatistics), verbose(verbose), fixedParameters(fixedParameters) {
}

////write a file with detailed stats for the optimization
//StaticParameterTuner::StaticParameterTuner(const std::string &tunerFileName, SGPP::base::OCLOperationConfiguration &fixedParameters, bool collectStatistics, bool verbose) :
//        collectStatistics(collectStatistics), verbose(verbose), fixedParameters(fixedParameters) {
//    this->readFromFile(tunerFileName);
//}

//void StaticParameterTuner::addFixedParameter(const std::string &name, const std::string &value, const ParameterType type) {
//    if (type == ParameterType::TEXT) {
//        this->kernelNode.addTextAttr(name, value);
//    } else if (type == ParameterType::ID) {
//        this->kernelNode.addIDAttr(name, value);
//    }
//}

void StaticParameterTuner::addParameter(const std::string &name, const std::vector<std::string> &valueRange) {
    this->tunableParameters.push_back(TunableParameter(name, valueRange));
}

SGPP::base::OCLOperationConfiguration StaticParameterTuner::tuneEverything(
SGPP::datadriven::LearnerScenario &scenario, const std::string &kernelName) {
    std::vector<std::string> platformsCopy = this->fixedParameters["PLATFORMS"].keys();
    for (const std::string &platformName : platformsCopy) {
        json::Node &platformNode = this->fixedParameters["PLATFORMS"][platformName];
        std::vector<std::string> devicesCopy = platformNode["DEVICES"].keys();

        //temporarily remove all other platforms
        std::vector<std::string> otherPlatformsName;
        std::vector<std::unique_ptr<json::Node>> otherPlatforms;
        for (const std::string &otherPlatformName : platformsCopy) {
            if (otherPlatformName.compare(platformName) == 0) {
                continue;
            }
            otherPlatformsName.push_back(otherPlatformName);
            auto otherPlatformNode = this->fixedParameters["PLATFORMS"][otherPlatformName].erase();
            otherPlatforms.push_back(std::move(otherPlatformNode));
        }

        for (const std::string &deviceName : devicesCopy) {
            json::Node &deviceNode = platformNode["DEVICES"][deviceName];

            std::cout << "tuning for device: " << deviceName << std::endl;

            //temporarily remove all other devices
            std::vector<std::string> otherDevicesName;
            std::vector<std::unique_ptr<json::Node>> otherDevices;
            for (const std::string &otherDeviceName : devicesCopy) {
                if (otherDeviceName.compare(deviceName) == 0) {
                    continue;
                }
                otherDevicesName.push_back(otherDeviceName);
                auto otherDeviceNode = platformNode["DEVICES"][otherDeviceName].erase();
                otherDevices.push_back(std::move(otherDeviceNode));
            }

            //limit the number of devices used to 1 for tuning
            bool addedDeviceLimit = false;
            size_t oldLimitValue;
            if (!deviceNode.contains("COUNT")) {
                deviceNode.addIDAttr("COUNT", 1ul);
                addedDeviceLimit = true;
            } else {
                oldLimitValue = deviceNode["COUNT"].getUInt();
            }

            //add an attribute for the kernel if none exists
            if (!deviceNode["KERNELS"].contains(kernelName)) {
                deviceNode["KERNELS"].addDictAttr(kernelName);
            }

            if (!deviceNode["KERNELS"][kernelName].contains("REUSE_SOURCE")) {
                deviceNode["KERNELS"][kernelName].addIDAttr("REUSE_SOURCE", false);
            }

            if (!deviceNode["KERNELS"][kernelName].contains("WRITE_SOURCE")) {
                deviceNode["KERNELS"][kernelName].addIDAttr("WRITE_SOURCE", false);
            }

//            if (useDoublePrecision) {
//                deviceNode["KERNELS"][kernelName].replaceTextAttr("INTERNAL_PRECISION", "double");
//            } else {
//                deviceNode["KERNELS"][kernelName].replaceTextAttr("INTERNAL_PRECISION", "float");
//            }

            this->tuneParameters(scenario, platformName, deviceName, kernelName);

            if (collectStatistics) {
                std::string safePlatformName = platformName;
                std::replace(safePlatformName.begin(), safePlatformName.end(), ' ', '_');
                std::string safeDeviceName = deviceName;
                std::replace(safeDeviceName.begin(), safeDeviceName.end(), ' ', '_');
                std::string statisticsFileName = "statistics_" + safePlatformName + "_" + safeDeviceName + "_"
                        + kernelName + ".csv";
                this->writeStatisticsToFile(statisticsFileName, platformName, deviceName, kernelName);
            }

            if (addedDeviceLimit) {
                deviceNode["COUNT"].erase();
            } else {
                deviceNode["COUNT"].setUInt(oldLimitValue);
            }

            //add the removed devices again for the next iteration
            for (size_t i = 0; i < otherDevices.size(); i++) {
                platformNode["DEVICES"].addAttribute(otherDevicesName[i], std::move(otherDevices[i]));
            }

        }

        //add the removed platforms again for the next iteration
        for (size_t i = 0; i < otherPlatforms.size(); i++) {
            this->fixedParameters["PLATFORMS"].addAttribute(otherPlatformsName[i], std::move(otherPlatforms[i]));
        }
    }

    std::cout << "after tuning: " << std::endl;
    std::cout << this->fixedParameters << std::endl;

    return this->fixedParameters;
}

void StaticParameterTuner::tuneParameters(SGPP::datadriven::LearnerScenario &scenario, const std::string &platformName,
        const std::string &deviceName, const std::string &kernelName) {

    if (collectStatistics) {
        this->statistics.clear();
    }

//    SGPP::base::OCLOperationConfiguration &currentParameters = fixedParameters;

    for (std::string &platformKey : fixedParameters["PLATFORMS"].keys()) {
        if (platformName.compare(platformKey) != 0) {
            throw;
        }
        for (std::string &deviceKey : fixedParameters["PLATFORMS"][platformKey]["DEVICES"].keys()) {
            if (deviceName.compare(deviceKey) != 0) {
                throw;
            }
            if (!fixedParameters["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"].contains(kernelName)) {
                throw;
            }
        }
    }

    json::Node &kernelNode = fixedParameters["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"][kernelName];

    //create initial parameter combination
    std::vector<size_t> valueIndices(tunableParameters.size());
    for (size_t i = 0; i < valueIndices.size(); i++) {
        valueIndices[i] = 0;
        TunableParameter &parameter = tunableParameters[i];
        kernelNode.replaceIDAttr(parameter.getName(), parameter.getValues()[0]);
    }

    std::cout << "-----------------------------------" << std::endl;
    for (std::string &key : kernelNode.keys()) {
        std::cout << "key: " << key << " value: " << kernelNode[key].get() << std::endl;
    }

    // evaluate initial parameter combination
    double shortestDuration = 0.0;//(scenario, fixedParameters, kernelName);
    if (collectStatistics) {
        this->statistics.push_back(std::make_pair(fixedParameters, shortestDuration));
    }

    std::unique_ptr<json::Node> bestParameters(kernelNode.clone());

    size_t parameterIndex = 0;
    while (parameterIndex < tunableParameters.size()) {

        // enumerate next parameter combination
        // can the current index be increased?
        TunableParameter &parameter = tunableParameters[parameterIndex];
        if (valueIndices[parameterIndex] + 1 < parameter.getValues().size()) {
            valueIndices[parameterIndex] += 1;
            kernelNode[parameter.getName()].set(parameter.getValues()[valueIndices[parameterIndex]]);
            //reset lower indices
            for (size_t i = 0; i < parameterIndex; i++) {
                valueIndices[i] = 0;
                TunableParameter &parameterForReset = tunableParameters[i];
                kernelNode[parameterForReset.getName()].set(parameterForReset.getValues()[0]);
            }
            parameterIndex = 0;

            std::cout << "-----------------------------------" << std::endl;
            for (std::string &key : kernelNode.keys()) {
                std::cout << "key: " << key << " value: " << kernelNode[key].get() << std::endl;
            }

            // evaluate current parameter combination
            double duration = 0.0;//(scenario, fixedParameters, kernelName);
            if (duration < shortestDuration) {
                std::cout << "new best combination! old: " << shortestDuration << " new: " << duration << std::endl;
                shortestDuration = duration;
                bestParameters = std::unique_ptr<json::Node>(kernelNode.clone());
                std::cout << *bestParameters << std::endl;
            }
            if (collectStatistics) {
                this->statistics.push_back(std::make_pair(fixedParameters, duration));
            }
        } else {
            parameterIndex += 1;
        }
    }

    std::cout << "overall best parameters:" << std::endl;
    std::cout << *bestParameters << std::endl;


//    this->fixedParameters["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"][kernelName] = *bestParameters;
    kernelNode = *bestParameters;

    std::cout << "written parameters:" << std::endl;
    std::cout << this->fixedParameters["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"][kernelName] << std::endl;
}

/*double StaticParameterTuner::0.0;//(SGPP::datadriven::LearnerScenario &scenario,
SGPP::base::OCLOperationConfiguration &currentParameters, const std::string &kernelName) {

    SGPP::datadriven::MetaLearner learner(scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
            scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(), scenario.getLambda(),
            verbose);

    SGPP::datadriven::OperationMultipleEvalType operationType;
    SGPP::datadriven::OperationMultipleEvalSubType operationSubType;

    //TODO: add better approach for auto-selecting the right kernel
    if (kernelName.compare("StreamingOCLMultiPlatform") == 0) {
        operationType = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
        operationSubType = SGPP::datadriven::OperationMultipleEvalSubType::OCLMP;
    } else if (kernelName.compare("StreamingModOCLFastMultiPlatform") == 0) {
        operationType = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
        operationSubType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;
    } else {
        throw;
    }

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(operationType, operationSubType,
            currentParameters);

    std::string fileName = scenario.getDatasetFileName();

    double duration = std::numeric_limits<double>::max();
    try {
        std::cout << "evaluating parameter combination" << std::endl;
        learner.learn(configuration, fileName);

        LearnerTiming timing = learner.getLearnerTiming();

        duration = timing.timeComplete_;
    } catch (SGPP::base::operation_exception &exception) {
        if (verbose) {
            std::cout << "invalid combination detected" << std::endl;
        }
    }

    if (verbose) {
        std::cout << "duration: " << duration << std::endl;
    }
    return duration;
}*/


/*void StaticParameterTuner::writeToFile(const std::string &fileName) {
 std::ofstream file(fileName);

 if (file.is_open()) {

 file << "fixed parameters" << std::endl;

 for (std::string &key : this->fixedParameters.keys()) {
 file << key << "=" << this->fixedParameters[key].get() << std::endl;
 }

 file << "tuned parameters" << std::endl;

 for (TunableParameter &parameter : this->tunableParameters) {

 file << parameter.getName() << "=";

 bool first = true;
 for (std::string &value : parameter.getValues()) {
 if (!first) {
 file << ",";
 } else {
 first = false;
 }
 file << value;
 }
 file << std::endl;
 }

 } else {
 throw;
 }

 file.close();
 }*/

/*void StaticParameterTuner::readFromFile(const std::string &fileName) {

 //reset instance if already initialized
 this->fixedParameters.clear();
 this->tunableParameters.clear();

 std::ifstream file(fileName);

 if (file.is_open()) {

 enum class ParserState {
 INITIAL, FIXED, TUNABLE
 };

 ParserState state = ParserState::INITIAL;

 std::string line;

 while (std::getline(file, line)) {

 std::vector<std::string> commentSplitted;
 boost::split(commentSplitted, line, boost::is_any_of("#"));
 std::string withoutComment = commentSplitted[0];
 boost::algorithm::trim(withoutComment);

 if (withoutComment.size() == 0) {
 continue;
 }

 if (state == ParserState::INITIAL) {
 if (withoutComment.compare("fixed parameters") == 0) {
 state = ParserState::FIXED;
 continue;
 } else {
 throw;
 }
 } else if (state == ParserState::FIXED) {
 if (withoutComment.compare("tuned parameters") == 0) {
 state = ParserState::TUNABLE;
 continue;
 }
 }

 std::vector<std::string> keyValueSplitted;
 boost::split(keyValueSplitted, withoutComment, boost::is_any_of("="));

 if (keyValueSplitted.size() != 2) {
 throw;
 }
 std::string key = keyValueSplitted[0];
 boost::algorithm::trim(key);
 std::string value = keyValueSplitted[1];
 boost::algorithm::trim(value);

 if (state == ParserState::FIXED) {
 this->addFixedParameter(key, value, ParameterType::ID);
 } else if (state == ParserState::TUNABLE) {
 std::vector<std::string> valuesSplitted;
 boost::split(valuesSplitted, value, boost::is_any_of(","));
 for (std::string &singleValue : valuesSplitted) {
 boost::algorithm::trim(singleValue);
 }
 this->addParameter(key, valuesSplitted, ParameterType::ID);
 } else {
 throw;
 }
 }

 } else {
 throw;
 }

 file.close();
 }*/

void StaticParameterTuner::writeStatisticsToFile(const std::string &statisticsFileName, const std::string &platformName,
        const std::string &deviceName, const std::string &kernelName) {
    if (!collectStatistics) {
        throw;
    }

    std::ofstream file(statisticsFileName);

    bool first = true;
    for (TunableParameter &columnParameter : this->tunableParameters) {
        if (!first) {
            file << ", ";
        } else {
            first = false;
        }
        file << columnParameter.getName();
    }
    if (this->tunableParameters.size() > 0) {
        file << ", ";
    }
    file << "duration" << std::endl;

    for (auto &parameterDurationPair : this->statistics) {
        SGPP::base::OCLOperationConfiguration&parameter = parameterDurationPair.first;
        json::Node &kernelNode = parameter["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"][kernelName];
        double duration = parameterDurationPair.second;

        first = true;
        for (TunableParameter &columnParameter : this->tunableParameters) {
            if (!first) {
                file << ", ";
            } else {
                first = false;
            }
            file << kernelNode[columnParameter.getName()].get();
        }
        if (this->tunableParameters.size() > 0) {
            file << ", ";
        }
        file << duration << std::endl;

    }
}

}
}

#endif
