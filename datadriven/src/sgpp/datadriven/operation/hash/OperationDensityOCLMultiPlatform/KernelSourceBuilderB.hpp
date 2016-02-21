
#pragma once

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename real_type>
class SourceBuilderB: public base::KernelSourceBuilderBase<real_type> {
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

	std::string generateSource(size_t dimensions, size_t datapoints) {
		if(kernelConfiguration.contains("REUSE_SOURCE")) {
			if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
				return this->reuseSource("DensityOCLMultiPlatform_mult.cl");
			}
		}

		std::stringstream sourceStream;

		if (this->floatType().compare("double") == 0) {
		sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
		}

		sourceStream<<"void kernel cscheme(global const int* starting_points,"<<std::endl
					<<"global const "<<this->floatType()<<"* data_points,global "<<this->floatType()<<"* C) {"<<std::endl
					<<"C[get_global_id(0)]=0.0;"<<std::endl
					<<"for(unsigned int ds=0;ds< " << datapoints << ";ds++)"<<std::endl
					<<"{"<<std::endl
					<<"private "<<this->floatType()<<" value=1;"<<std::endl
					<<"for(private int d=0;d< " << dimensions << ";d++)"<<std::endl
					<<"{"<<std::endl
					<<"private "<<this->floatType()<<" wert=1.0;"<<std::endl
					<<"for(private int z=1;"<<std::endl
					<<"	z<=starting_points[(get_global_id(0))*2* " << dimensions << "+2*d+1];z++)"<<std::endl
					<<"		wert*=2;"<<std::endl
					<<"wert*=data_points[ds* " << dimensions << "+d];"<<std::endl
					<<"	wert-=starting_points[(get_global_id(0))*2* " << dimensions << "+2*d];"<<std::endl
					<<"	if(wert<0)"<<std::endl
					<<"		wert*=-1;"<<std::endl
					<<"	wert=1-wert;"<<std::endl
					<<"	if(wert<0)"<<std::endl
					<<"		wert=0;"<<std::endl
					<<"	value*=wert;"<<std::endl
					<<"}"<<std::endl
					<<"C[get_global_id(0)]+=value/ " << datapoints << ";"<<std::endl
					<<"}"<<std::endl
					<<"}"<<std::endl;
		if(kernelConfiguration.contains("WRITE_SOURCE")) {
			if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
				this->writeSource("DensityOCLMultiPlatform_mult.cl", sourceStream.str());
			}
		}


		return sourceStream.str();
	}

}
	;

}
}
}
