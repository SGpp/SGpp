/*
 * test_datadrivenCommon.cpp
 *
 *  Created on: Jul 30, 2015
 *      Author: pfandedd
 */

#include "test_datadrivenCommon.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include <zlib.h>

using namespace SGPP::base;

std::string uncompressFile(std::string fileName) {

  gzFile inFileZ = gzopen(fileName.c_str(), "rb");

  if (inFileZ == NULL) {
    std::cout << "Error: Failed to gzopen file " << fileName << std::endl;
    exit(0);
  }

  unsigned char unzipBuffer[8192];
  unsigned int unzippedBytes;
  std::vector<unsigned char> unzippedData;

  while (true) {
    unzippedBytes = gzread(inFileZ, unzipBuffer, 8192);

    if (unzippedBytes > 0) {
      for (size_t i = 0; i < unzippedBytes; i++) {
        unzippedData.push_back(unzipBuffer[i]);
      }
    } else {
      break;
    }
  }

  gzclose(inFileZ);

  std::stringstream convert;

  for (size_t i = 0; i < unzippedData.size(); i++) {
    convert << unzippedData[i];
  }

  return convert.str();
}

DataMatrix* readReferenceMatrix(GridStorage* storage, std::string fileName) {

  std::string content = uncompressFile(fileName);

  std::stringstream contentStream;
  contentStream << content;
  std::string line;

  DataMatrix* m = new DataMatrix(0, storage->size());

  size_t currentRow = 0;

  while (!contentStream.eof()) {

    std::getline(contentStream, line);

    // for lines that only contain a newline
    if (line.size() == 0) {
      break;
    }

    m->appendRow();

    size_t curPos = 0;
    size_t curFind = 0;
    std::string curValue;
    float_t floatValue;

    for (size_t i = 0; i < storage->size(); i++) {
      curFind = line.find_first_of(" \t", curPos);
      curValue = line.substr(curPos, curFind - curPos);

      floatValue = boost::lexical_cast<float_t>(curValue);
      m->set(currentRow, i, floatValue);
      curPos = curFind + 1;
    }

    currentRow += 1;
  }

  return m;
}
