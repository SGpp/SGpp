/*
 * StreamingOCLParameters.hpp
 *
 *  Created on: Mar 18, 2015
 *      Author: pfandedd
 */


#ifndef STREAMING_OCL_DEVICE_TYPE
//options: CL_DEVICE_TYPE_CPU, CL_DEVICE_TYPE_GPU, ...
//this define requires the include of <CL/cl.h>
#define STREAMING_OCL_DEVICE_TYPE CL_DEVICE_TYPE_CPU
#endif

#ifndef STREAMING_OCL_ENABLE_OPTIMIZATIONS
#define STREAMING_OCL_ENABLE_OPTIMIZATIONS true
#endif

#ifndef STREAMING_OCL_USE_LOCAL_MEMORY
#define STREAMING_OCL_USE_LOCAL_MEMORY true
#endif

#ifndef STREAMING_OCL_LOCAL_SIZE
#define STREAMING_OCL_LOCAL_SIZE 64
#endif
