
#pragma once

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename real_type>
class SourceBuilderCreateGraph: public base::KernelSourceBuilderBase<real_type> {
private:

	std::shared_ptr<base::OCLDevice> device;

	json::Node &kernelConfiguration;

	size_t dims;

	size_t localWorkgroupSize;
	bool useLocalMemory;
	size_t dataBlockSize;
	size_t transGridBlockSize;
	uint64_t maxDimUnroll;

	std::string getData(std::string dim, size_t dataBlockingIndex) {
		std::stringstream output;
		if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("array") == 0) {
			output << "data_" << dataBlockingIndex << "[" << dim << "]";
		} else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
			output << "data_" << dataBlockingIndex << "_" << dim;
		} else if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("pointer") == 0) {
			output << "ptrData[(" << dataBlockSize << " * globalIdx) + (resultSize * " << dim << ") + "
				   << dataBlockingIndex << "]";
		} else {
			throw new base::operation_exception("OCL error: Illegal value for parameter \"KERNEL_STORE_DATA\"\n");
		}
		return output.str();
	}
public:

	SourceBuilderCreateGraph(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration, size_t dims) :
		device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
	}

	std::string generateSource(int dimensions, int k, int dataSize) {
		if(kernelConfiguration.contains("REUSE_SOURCE")) {
			if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
				return this->reuseSource("DensityOCLMultiPlatform_mult.cl");
			}
		}

		std::stringstream sourceStream;

		if (this->floatType().compare("double") == 0) {
			sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
		}

		sourceStream<<""<<std::endl
					<<"void kernel connectNeighbors(global "<<this->floatType()<<" *data,"<<std::endl
					<<"							    global int *neighbors)"<<std::endl
					<<"{"<<std::endl
					<<"	  size_t global_index = get_global_id(0);"<<std::endl
					<<"	  size_t index = get_global_id(0);"<<std::endl
					<<"	  for(int t=0;t< " << k << " ;t++)"<<std::endl
					<<"		  neighbors[( " << k << " *index)+t] = t;"<<std::endl
					<<"	  int size = 0;"<<std::endl
					<<"	  for (unsigned int i = 0; i <  " << dataSize << "; i++) {"<<std::endl
					<<"		 if (i != global_index) {"<<std::endl
					<<"			 //get distance to current point"<<std::endl
					<<"			"<<this->floatType()<<" dist = 0.0;"<<std::endl
					<<"			 for (unsigned int j = 0; j <  " << dimensions << " ; j++) {"<<std::endl
					<<"			 dist += (data[global_index* " << dimensions << "  + j] - data[j + i* " << dimensions << " ])"<<std::endl
					<<"				   * (data[j + global_index* " << dimensions << " ] - data[j + i* " << dimensions << " ]);"<<std::endl
					<<"		 }"<<std::endl
					<<"			"<<this->floatType()<<" max=0.0;"<<std::endl
					<<"			int maxindex=-1;"<<std::endl
					<<"			for(int t=0;t< " << k << " ;t++)"<<std::endl
					<<"			{"<<std::endl
					<<"				int currentneighbor=neighbors[ " << k << " *index+t];"<<std::endl
					<<"				"<<this->floatType()<<" currentdist=0.0;"<<std::endl
					<<"				for (unsigned int j = 0; j <  " << dimensions << " ; j++) {"<<std::endl
					<<"				   currentdist += (data[global_index* " << dimensions << "  + j] - data[j + currentneighbor* " << dimensions << " ])"<<std::endl
					<<"								* (data[j + global_index* " << dimensions << " ] - data[j + currentneighbor* " << dimensions << " ]);"<<std::endl
					<<"				}"<<std::endl
					<<"				if(max<currentdist)"<<std::endl
					<<"				{"<<std::endl
					<<"				   max=currentdist;"<<std::endl
					<<"				   maxindex=t;"<<std::endl
					<<"				}"<<std::endl
					<<"			 }"<<std::endl
					<<"			 if(dist<max)"<<std::endl
					<<"			 {"<<std::endl
					<<"				neighbors[( " << k << " *index)+maxindex] = i;"<<std::endl
					<<"			 }"<<std::endl
					<<"		  }"<<std::endl
					<<"	   }"<<std::endl
					<<"	}"<<std::endl;
		if(kernelConfiguration.contains("WRITE_SOURCE")) {
			if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
				this->writeSource("DensityOCLMultiPlatform_graph.cl", sourceStream.str());
			}
		}


		return sourceStream.str();
	}

}
	;

}
}
}
