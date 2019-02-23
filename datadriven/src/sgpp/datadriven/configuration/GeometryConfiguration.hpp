/*
 * GeometryConfiguration.hpp
 *
 *  Created on: Jan 15, 2019
 *      Author: jamal
 */

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_CONFIGURATION_GEOMETRYCONFIGURATION_HPP_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_CONFIGURATION_GEOMETRYCONFIGURATION_HPP_

#include <sgpp/globaldef.hpp>
#include <vector>
#include <string>


namespace sgpp {
namespace datadriven {

/*
 * Struct that stores information to geometry aware sparse grids
 */
struct GeometryConfiguration{
  /*
   * Stencil for geometric relation
   */
  std::string stencil;

  /*
   * resolution of image/video e.g 28x28
   */
  std::vector<int64_t> dim;
};


}  // namespace datadriven
}  // namespace sgpp




#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_CONFIGURATION_GEOMETRYCONFIGURATION_HPP_ */
