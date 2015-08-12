// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

Dataset ARFFTools::readARFF(const std::string& filename) {
    std::string line;
    std::ifstream myfile(filename.c_str());
    size_t numberInstances;
    size_t dimension;
    bool dataReached = false;
    size_t instanceNo = 0;

    readARFFSize(filename, numberInstances, dimension);
    Dataset dataset(numberInstances, dimension);

    while (!myfile.eof()) {
        std::getline(myfile, line);
        std::transform(line.begin(), line.end(), line.begin(), toupper);

        if (dataReached && !line.empty()) {
            writeNewClass(line, *dataset.getClasses(), instanceNo);
            writeNewTrainingDataEntry(line, *dataset.getTrainingData(), instanceNo);
            instanceNo++;
        }

        if (line.find("@DATA", 0) != line.npos) {
            dataReached = true;
        }
    }

    myfile.close();

    return dataset;
}

void ARFFTools::readARFFSize(const std::string& filename, size_t& numberInstances, size_t& dimension) {
    std::string line;
    std::ifstream myfile(filename.c_str());
    dimension = 0;
    numberInstances = 0;

    if (!myfile.is_open()) {
        std::string msg = "Unable to open file: " + filename;
        throw new SGPP::base::file_exception(msg.c_str());
    }

    while (!myfile.eof()) {
        std::getline(myfile, line);
        std::transform(line.begin(), line.end(), line.begin(), toupper);

        if (line.find("@ATTRIBUTE class", 0) != line.npos) {
            ;
        } else if (line.find("@ATTRIBUTE CLASS", 0) != line.npos) {
            ;
        } else if (line.find("@ATTRIBUTE", 0) != line.npos) {
            dimension++;
        } else if (line.find("@DATA", 0) != line.npos) {
            numberInstances = 0;
        } else if (!line.empty()) {
            numberInstances++;
        }
    }

    myfile.close();
}

void ARFFTools::readARFFSizeFromString(const std::string& content, size_t& numberInstances, size_t& dimension) {
    std::string line;
    std::istringstream contentStream(content);

    dimension = 0;
    numberInstances = 0;

    while (!contentStream.eof()) {
        std::getline(contentStream, line);
        std::transform(line.begin(), line.end(), line.begin(), toupper);

        if (line.find("@ATTRIBUTE class", 0) != line.npos) {
            ;
        } else if (line.find("@ATTRIBUTE CLASS", 0) != line.npos) {
            ;
        } else if (line.find("@ATTRIBUTE", 0) != line.npos) {
            dimension++;
        } else if (line.find("@DATA", 0) != line.npos) {
            numberInstances = 0;
        } else if (!line.empty()) {
            numberInstances++;
        }
    }

}

Dataset ARFFTools::readARFFFromString(const std::string& content) {
    std::string line;
    std::stringstream contentStream;
    contentStream << content;
    size_t numberInstances;
    size_t dimension;
    bool dataReached = false;
    size_t instanceNo = 0;

    ARFFTools::readARFFSizeFromString(content, numberInstances, dimension);
    Dataset dataset(numberInstances, dimension);

    while (!contentStream.eof()) {
        std::getline(contentStream, line);
        std::transform(line.begin(), line.end(), line.begin(), toupper);

        if (dataReached && !line.empty()) {
            writeNewClass(line, *dataset.getClasses(), instanceNo);
            writeNewTrainingDataEntry(line, *dataset.getTrainingData(), instanceNo);
            instanceNo++;
        }

        if (line.find("@DATA", 0) != line.npos) {
            dataReached = true;
        }
    }

    return dataset;
}

void ARFFTools::writeNewTrainingDataEntry(const std::string& arffLine, SGPP::base::DataMatrix& destination,
        size_t instanceNo) {
    size_t cur_pos = 0;
    size_t cur_find = 0;
    size_t dim = destination.getNcols();
    std::string cur_value;
    float_t dbl_cur_value;

    for (size_t i = 0; i < dim; i++) {
        cur_find = arffLine.find(",", cur_pos);
        cur_value = arffLine.substr(cur_pos, cur_find - cur_pos);
        dbl_cur_value = atof(cur_value.c_str());
        destination.set(instanceNo, i, dbl_cur_value);
        cur_pos = cur_find + 1;
    }
}

void ARFFTools::writeNewClass(const std::string& arffLine, SGPP::base::DataVector& destination, size_t instanceNo) {
    size_t cur_pos = arffLine.find_last_of(",");
    std::string cur_value = arffLine.substr(cur_pos + 1);
    float_t dbl_cur_value = atof(cur_value.c_str());
    destination.set(instanceNo, dbl_cur_value);
}

//void ARFFTools::writeAlpha(std::string tfilename, SGPP::base::DataVector& source)
//{
//
//}

//void ARFFTools::readAlpha(std::string tfilename, SGPP::base::DataVector& destination)
//{
//
//}

}
}
