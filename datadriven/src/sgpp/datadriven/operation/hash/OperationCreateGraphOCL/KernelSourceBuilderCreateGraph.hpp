
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

	std::string getData(size_t dim, size_t dataBlockingIndex) {
		std::stringstream dimStringStream;
		dimStringStream << dim;
		std::string dimString = dimStringStream.str();
		return this->getData(dimString, dataBlockingIndex);
	}

	std::string unrolledBasisFunctionEvalulation(size_t dims, size_t startDim, size_t endDim,
												 std::string unrollVariable) {
		std::stringstream output;

		for (size_t d = startDim; d < endDim; d++) {

			std::stringstream dimElement;
			dimElement << "(";
			if (!unrollVariable.compare("") == 0) {
				dimElement << unrollVariable << " + ";
			}
			dimElement << d;
			dimElement << ")";
			std::string pointerAccess = dimElement.str();

			std::string dString;
			if (kernelConfiguration["KERNEL_STORE_DATA"].get().compare("register") == 0) {
				std::stringstream stream;
				stream << (d);
				dString = stream.str();
			} else {
				dString = pointerAccess;
			}

			std::stringstream levelAccessStream;
			std::stringstream indexAccessStream;
			if (useLocalMemory) {
				levelAccessStream << "locLevel[dimLevelIndex]";
				indexAccessStream << "locIndex[dimLevelIndex]";
			} else {
				levelAccessStream << "ptrLevel[dimLevelIndex]";
				indexAccessStream << "ptrIndex[dimLevelIndex]";
			}
			std::string levelAccess = levelAccessStream.str();
			std::string indexAccess = indexAccessStream.str();

			output << this->indent3 << "dimLevelIndex = " << "(k * " << dims << ") + " << pointerAccess << ";"
				   << std::endl;

			for (size_t i = 0; i < dataBlockSize; i++) {
				output << this->indent3 << "curSupport_" << i << " *= fmax(1.0" << this->constSuffix() << " - fabs((";
				output << levelAccess << " * " << getData(dString, i) << ") - " << indexAccess << "), 0.0"
					   << this->constSuffix() << ");" << std::endl;
			}
		}
		return output.str();
	}

public:

	SourceBuilderB(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration, size_t dims) :
		device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
	}

	std::string generateSource() {
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
					<<"__kernel void connectNeighbors(__global "<<this->floatType()<<" *data, const int dimensions,"<<std::endl
					<<"							   const int k, const int dataSize, __global int *neighbors)"<<std::endl
					<<"{"<<std::endl
					<<"	size_t global_index = get_global_id(0);"<<std::endl
					<<"	size_t index = get_global_id(0);"<<std::endl
					<<"	int size = 0;"<<std::endl
					<<"	for (unsigned int i = 0; i < dataSize/dimensions; i++) {"<<std::endl
					<<"	if (i != global_index) {"<<std::endl
					<<"		//get distance to current point"<<std::endl
					<<"	"<<this->floatType()<<" dist = 0.0;"<<std::endl
					<<"	for (unsigned int j = 0; j < dimensions; j++) {"<<std::endl
					<<"		dist += (data[global_index*dimensions + j] - data[j + i*dimensions])"<<std::endl
					<<"			* (data[j + global_index*dimensions] - data[j + i*dimensions]);"<<std::endl
					<<"	}"<<std::endl
					<<""<<std::endl
					<<"			//Not enough neighbors yet"<<std::endl
					<<"			if (size < k)"<<std::endl
					<<"		{"<<std::endl
					<<"			//Just add current point as a neighbori"<<std::endl
					<<"		neighbors[(k*index)+size] = i;"<<std::endl
					<<"		size++;"<<std::endl
					<<"	}"<<std::endl
					<<"	else"<<std::endl
					<<"		{"<<std::endl
					<<"		"<<this->floatType()<<" max=0.0;"<<std::endl
					<<"	int maxindex=-1;"<<std::endl
					<<"	for(int t=0;t<k;t++)"<<std::endl
					<<"		{"<<std::endl
					<<"			int currentneighbor=neighbors[k*index+t];"<<std::endl
					<<"		"<<this->floatType()<<" currentdist=0.0;"<<std::endl
					<<"		for (unsigned int j = 0; j < dimensions; j++) {"<<std::endl
					<<"			currentdist += (data[global_index*dimensions + j] - data[j + currentneighbor*dimensions])"<<std::endl
					<<"				* (data[j + global_index*dimensions] - data[j + currentneighbor*dimensions]);"<<std::endl
					<<""<<std::endl
					<<"if(max<currentdist)"<<std::endl
					<<"	{"<<std::endl
					<<"		max=currentdist;"<<std::endl
					<<"		maxindex=t;"<<std::endl
					<<"	}"<<std::endl
					<<"	}"<<std::endl
					<<"	if(dist<max)"<<std::endl
					<<"	{"<<std::endl
					<<"		neighbors[(k*index)+maxindex] = i;"<<std::endl
					<<"	}"<<std::endl
					<<"	}"<<std::endl
					<<"	}"<<std::endl
					<<"	}"<<std::endl
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
