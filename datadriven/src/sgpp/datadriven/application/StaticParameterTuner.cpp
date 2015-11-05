#if USE_OCL == 1

#include <limits>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/datadriven/application/StaticParameterTuner.hpp>
#include "sgpp/datadriven/application/MetaLearner.hpp"

namespace SGPP {
namespace datadriven {

StaticParameterTuner::StaticParameterTuner(bool collectStatistics, bool verbose) :
        collectStatistics(collectStatistics), verbose(verbose) {

}

//write a file with detailed stats for the optimization
StaticParameterTuner::StaticParameterTuner(std::string tunerFileName, bool collectStatistics, bool verbose) :
        collectStatistics(collectStatistics), verbose(verbose) {
    this->readFromFile(tunerFileName);
}

void StaticParameterTuner::addFixedParameter(std::string name, std::string value) {
    this->fixedParameters.set(name, value);
}

void StaticParameterTuner::addParameter(std::string name, std::vector<std::string> valueRange) {
    this->tunableParameters.push_back(TunableParameter(name, valueRange));
}

SGPP::base::OCLConfigurationParameters StaticParameterTuner::tuneParameters(
SGPP::datadriven::LearnerScenario scenario) {

    if (collectStatistics) {
        this->statistics.clear();
    }

    SGPP::base::OCLConfigurationParameters currentParameters = fixedParameters;

    //create initial parameter combination
    std::vector<size_t> valueIndices(tunableParameters.size());
    for (size_t i = 0; i < valueIndices.size(); i++) {
        valueIndices[i] = 0;
        TunableParameter &parameter = tunableParameters[i];
        currentParameters.set(parameter.getName(), parameter.getValues()[0]);
    }

    // evaluate initial parameter combination
    double shortestDuration = evaluateSetup(scenario, currentParameters);
    SGPP::base::OCLConfigurationParameters bestParameters = currentParameters;

    size_t parameterIndex = 0;
    while (parameterIndex < tunableParameters.size()) {

        // enumerate next parameter combination
        // can the current index be increased?
        TunableParameter &parameter = tunableParameters[parameterIndex];
        if (valueIndices[parameterIndex] + 1 < parameter.getValues().size()) {
            valueIndices[parameterIndex] += 1;
            currentParameters.set(parameter.getName(), parameter.getValues()[valueIndices[parameterIndex]]);
            //reset lower indices
            for (size_t i = 0; i < parameterIndex; i++) {
                valueIndices[i] = 0;
                TunableParameter &parameterForReset = tunableParameters[i];
                currentParameters.set(parameterForReset.getName(), parameterForReset.getValues()[0]);
            }
            parameterIndex = 0;

            // evaluate current parameter combination
            double duration = evaluateSetup(scenario, currentParameters);
            if (duration < shortestDuration) {
                shortestDuration = duration;
                bestParameters = currentParameters;
            }
            if (collectStatistics) {
                this->statistics.push_back(std::make_pair(currentParameters, duration));
            }
        } else {
            parameterIndex += 1;
        }
    }

    return bestParameters;
}

double StaticParameterTuner::evaluateSetup(SGPP::datadriven::LearnerScenario scenario,
SGPP::base::OCLConfigurationParameters currentParameters) {

    SGPP::datadriven::MetaLearner learner(scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
            scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(), scenario.getLambda(),
            verbose);

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, currentParameters);

    std::string fileName = scenario.getDatasetFileName();

    double duration = std::numeric_limits<double>::max();
    try {
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
}

void StaticParameterTuner::writeToFile(std::string fileName) {
    std::ofstream file(fileName);

    if (file.is_open()) {

        file << "fixed parameters" << std::endl;

        for (std::string key : this->fixedParameters.getKeys()) {
            file << key << "=" << this->fixedParameters.get(key) << std::endl;
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
}

void StaticParameterTuner::readFromFile(std::string fileName) {

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
                this->addFixedParameter(key, value);
            } else if (state == ParserState::TUNABLE) {
                std::vector<std::string> valuesSplitted;
                boost::split(valuesSplitted, value, boost::is_any_of(","));
                for (std::string &singleValue : valuesSplitted) {
                    boost::algorithm::trim(singleValue);
                }
                this->addParameter(key, valuesSplitted);
            } else {
                throw;
            }
        }

    } else {
        throw;
    }

    file.close();
}

void StaticParameterTuner::writeStatisticsToFile(std::string statisticsFileName) {
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

    for (auto parameterDurationPair : this->statistics) {
        SGPP::base::OCLConfigurationParameters &parameter = parameterDurationPair.first;
        double duration = parameterDurationPair.second;

        first = true;
        for (TunableParameter &columnParameter : this->tunableParameters) {
            if (!first) {
                file << ", ";
            } else {
                first = false;
            }
            file << parameter.get(columnParameter.getName());
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
