/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * ARFFWrapper.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun, Michael Lettrich
 */

#include <sgpp/datadriven/datamining/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <string>

namespace sgpp {

namespace datadriven {

ArffFileSampleProvider::ArffFileSampleProvider() {}

ArffFileSampleProvider::~ArffFileSampleProvider() {}

size_t ArffFileSampleProvider::getDim() { return 0; }

size_t ArffFileSampleProvider::getDatasetSize() { return 0; }

void ArffFileSampleProvider::readFile(std::string fileName) {}

std::unique_ptr<Dataset> ArffFileSampleProvider::getNextSamples(int howMany) { return nullptr; }

std::unique_ptr<Dataset> ArffFileSampleProvider::getAllSamples() { return nullptr; }

size_t ArffFileSampleProvider::getNumClasses() { return 0; }

void ArffFileSampleProvider::readString(std::string& input) {}

} /* namespace datadriven */
} /* namespace sgpp */
