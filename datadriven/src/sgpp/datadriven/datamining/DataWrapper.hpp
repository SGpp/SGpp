/*
 * DataWrapper.h
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef DATAWRAPPER_H_
#define DATAWRAPPER_H_

#include <string>

#include <sgpp/datadriven/datamining/SampleProvider.hpp>
#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class DataWrapper : public SampleProvider {
 public:
  DataWrapper(datadriven::DataMiningConfiguration& config) : SampleProvider(config) {
    try {
      filename = config["filename"].get();
    } catch (json::json_exception& e) {
      std::cout << e.what() << std::endl;
    }
  }
  virtual ~DataWrapper(){};

 protected:
  std::string filename;
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* DATAWRAPPER_H_ */
