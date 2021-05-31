// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/tools/StringTokenizer.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceFileTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/FileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataShufflingFunctorFactory.hpp>

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

using sgpp::base::data_exception;

DataSourceBuilder& DataSourceBuilder::withFileType(DataSourceFileType fileType) {
  config.fileType_ = fileType;
  return *this;
}

DataSourceBuilder& DataSourceBuilder::inBatches(size_t howMany) {
  config.numBatches_ = howMany;
  return *this;
}

DataSourceBuilder& DataSourceBuilder::withBatchSize(size_t batchSize) {
  config.batchSize_ = batchSize;
  return *this;
}

DataSourceBuilder& DataSourceBuilder::withCompression(bool isCompressed) {
  config.isCompressed_ = isCompressed;

#ifndef ZLIB
  if (isCompressed) {
    throw sgpp::base::application_exception{
        "sgpp has been built without zlib support. Reading compressed files is not possible"};
  }
#endif

  return *this;
}

DataSourceBuilder& DataSourceBuilder::withPath(const std::string& filePath) {
  config.filePath_ = filePath;
  if (config.fileType_ == DataSourceFileType::NONE) {
    grabTypeInfoFromFilePath();
  }
  return *this;
}

DataSourceSplitting* DataSourceBuilder::splittingAssemble() const {
  // Create a shuffling functor
  DataShufflingFunctorFactory shufflingFunctorFactory;
  DataShufflingFunctor* shuffling = shufflingFunctorFactory.buildDataShufflingFunctor(config);

  SampleProvider* sampleProvider = nullptr;

  if (config.fileType_ == DataSourceFileType::ARFF) {
    sampleProvider = new ArffFileSampleProvider(shuffling);
  } else if (config.fileType_ == DataSourceFileType::CSV) {
    sampleProvider = new CSVFileSampleProvider(shuffling);
  } else {
    throw data_exception("DataSourceBuilder::splittingAssemble() unknown file type");
  }

  if (config.isCompressed_) {
#ifndef ZLIB
    throw sgpp::base::application_exception{
        "sgpp has been built without zlib support. Reading compressed files is not possible"};
#else
    sampleProvider = new GzipFileSampleDecorator(static_cast<FileSampleProvider*>(sampleProvider));
#endif
  }

  return new DataSourceSplitting(config, sampleProvider);
}

DataSourceSplitting* DataSourceBuilder::splittingFromConfig(const DataSourceConfig& config) {
  this->config = config;

  if (config.fileType_ == DataSourceFileType::NONE) {
    grabTypeInfoFromFilePath();
  }
  return splittingAssemble();
}

DataSourceCrossValidation* DataSourceBuilder::crossValidationAssemble() const {
  // Create a shuffling functor
  DataShufflingFunctorFactory shufflingFunctorFactory;
  DataShufflingFunctor* shuffling = shufflingFunctorFactory.buildDataShufflingFunctor(config);
  DataShufflingFunctorCrossValidation* crossValidationShuffling =
      new DataShufflingFunctorCrossValidation(crossValidationConfig, shuffling);

  SampleProvider* sampleProvider = nullptr;

  if (config.fileType_ == DataSourceFileType::ARFF) {
    sampleProvider = new ArffFileSampleProvider(crossValidationShuffling);
  } else if (config.fileType_ == DataSourceFileType::CSV) {
    sampleProvider = new CSVFileSampleProvider(crossValidationShuffling);
  } else {
    throw data_exception("DataSourceBuilder::crossValidationAssemble() unknown file type");
  }

  if (config.isCompressed_) {
#ifndef ZLIB
    throw sgpp::base::application_exception{
        "sgpp has been built without zlib support. Reading compressed files is not possible"};
#else
    sampleProvider = new GzipFileSampleDecorator(static_cast<FileSampleProvider*>(sampleProvider));
#endif
  }
  auto t = new DataSourceCrossValidation(config, crossValidationConfig, crossValidationShuffling,
                                         sampleProvider);
  return t;
}

DataSourceCrossValidation* DataSourceBuilder::crossValidationFromConfig(
    const DataSourceConfig& config, const CrossvalidationConfiguration& crossValidationConfig) {
  this->config = config;
  this->crossValidationConfig = crossValidationConfig;

  if (config.fileType_ == DataSourceFileType::NONE) {
    grabTypeInfoFromFilePath();
  }
  return crossValidationAssemble();
}

void DataSourceBuilder::grabTypeInfoFromFilePath() {
  // tokenize string
  std::vector<std::string> tokens;
  // split the string
  sgpp::base::StringTokenizer::tokenize(config.filePath_, ".", tokens);
  // convert to lower case
  for (auto t : tokens) {
    // TODO(Michael Lettrich): test if this works with umlauts
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
  }
  // check if there is gz compression
  if (tokens.back() == "gz") {
    withCompression(true);
  }

  // check if we can find file type
  auto type = DataSourceFileType::NONE;
  for (auto t : tokens) {
    try {
      type = DataSourceFileTypeParser::parse(t);
    } catch (data_exception&) {
      // wasn't found
      withFileType(DataSourceFileType::NONE);
    }
    if (type != DataSourceFileType::NONE) {
      withFileType(type);
    }
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
