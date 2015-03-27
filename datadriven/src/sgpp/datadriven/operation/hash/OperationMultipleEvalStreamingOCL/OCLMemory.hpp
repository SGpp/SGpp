/*
 * OCLMemory.hpp
 *
 *  Created on: Mar 27, 2015
 *      Author: pfandedd
 */

namespace SGPP {
namespace datadriven {

struct OCLMemory {
  cl_mem *bufferList;
  size_t sizeofType;
  size_t elements;

  bool isMappedMemory;
  cl_mem hostBuffer;
  void *mappedHostBuffer;
};

}
}


