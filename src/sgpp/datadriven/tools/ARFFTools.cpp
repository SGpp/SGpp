/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "datadriven/tools/ARFFTools.hpp"
#include "base/exception/file_exception.hpp"
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

namespace sg {
  namespace datadriven {

    ARFFTools::ARFFTools() {
    }

    ARFFTools::~ARFFTools() {
    }

    size_t ARFFTools::getDimension(std::string tfilename) {
      std::string line;
      std::ifstream myfile (tfilename.c_str());
      size_t numFound = 0;

      if (myfile.is_open()) {
        while (!myfile.eof() ) {
          getline (myfile, line);
          std::transform(line.begin(), line.end(), line.begin(), toupper);
          if (line.find("@ATTRIBUTE", 0) != line.npos) {
            numFound++;
          }
        }

        myfile.close();
      } else {
        std::string msg = "Unable to open file: " + tfilename;
        throw new sg::base::file_exception(msg.c_str());
      }

      // the class is not regarded when getting the dimension
      numFound--;

      return numFound;
    }

    size_t ARFFTools::getNumberInstances(std::string tfilename) {
      std::string line;
      std::ifstream myfile (tfilename.c_str());
      size_t numInst = 0;

      if (myfile.is_open()) {
        getline (myfile, line);

        while (!myfile.eof() ) {
          std::transform(line.begin(), line.end(), line.begin(), toupper);
          if (line.find("@DATA", 0) != line.npos) {
            numInst = 0;
          } else {
            numInst++;
          }

          getline (myfile, line);
        }

        myfile.close();
      } else {
        std::string msg = "Unable to open file: " + tfilename;
        throw new sg::base::file_exception(msg.c_str());
      }

      return numInst;
    }

    void ARFFTools::readTrainingData(std::string tfilename, sg::base::DataMatrix& destination) {
      std::string line;
      std::ifstream myfile (tfilename.c_str());
      bool data = false;
      size_t instanceNo = 0;

      if (myfile.is_open()) {
        getline (myfile, line);

        while (!myfile.eof() ) {
          std::transform(line.begin(), line.end(), line.begin(), toupper);
          if (data == true) {
            writeNewElement(line, destination, instanceNo);
            instanceNo++;
          }

          if (line.find("@DATA", 0) != line.npos) {
            data = true;
          }

          getline (myfile, line);
        }

        myfile.close();
      } else {
        std::string msg = "Unable to open file: " + tfilename;
        throw new sg::base::file_exception(msg.c_str());
      }
    }

    void ARFFTools::readClasses(std::string tfilename, sg::base::DataVector& destination) {
      std::string line;
      std::ifstream myfile (tfilename.c_str());
      bool data = false;
      size_t instanceNo = 0;

      if (myfile.is_open()) {
        getline (myfile, line);

        while (!myfile.eof() ) {
          std::transform(line.begin(), line.end(), line.begin(), toupper);
          if (data == true) {
            writeNewClass(line, destination, instanceNo);
            instanceNo++;
          }

          if (line.find("@DATA", 0) != line.npos) {
            data = true;
          }

          getline (myfile, line);
        }

        myfile.close();
      } else {
        std::string msg = "Unable to open file: " + tfilename;
        throw new sg::base::file_exception(msg.c_str());
      }
    }

    void ARFFTools::writeNewElement(std::string& instance, sg::base::DataMatrix& destination, size_t instanceNo) {
      size_t cur_pos = 0;
      size_t cur_find = 0;
      size_t dim = destination.getNcols();
      std::string cur_value;
      double dbl_cur_value;

      for (size_t i = 0; i < dim; i++) {
        cur_find = instance.find(",", cur_pos);
        cur_value = instance.substr(cur_pos, cur_find - cur_pos);
        dbl_cur_value = atof(cur_value.c_str());
        destination.set(instanceNo , i, dbl_cur_value);
        cur_pos = cur_find + 1;
      }
    }

    void ARFFTools::writeNewClass(std::string& instance, sg::base::DataVector& destination, size_t instanceNo) {
      size_t cur_pos = instance.find_last_of(",");
      std::string cur_value = instance.substr(cur_pos + 1);
      double dbl_cur_value = atof(cur_value.c_str());
      destination.set(instanceNo, dbl_cur_value);
    }

    //void ARFFTools::writeAlpha(std::string tfilename, sg::base::DataVector& source)
    //{
    //
    //}

    //void ARFFTools::readAlpha(std::string tfilename, sg::base::DataVector& destination)
    //{
    //
    //}

  }
}
