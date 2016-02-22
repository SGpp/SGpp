
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

	std::string init_k_registers(size_t k, size_t dimensions) {
		std::stringstream output;
		for(size_t i = 0; i < k; i++) {
			output << this->indent[0] << "int k_register" << i << " = " << i << "; " << std::endl;
		}
		for(size_t i = 0; i < k; i++) {
			output << this->indent[0] << this->floatType() << " dist_k" << i << " = 4.0;"<<std::endl;
			/*output << "for (unsigned int j = 0; j <  " << dimensions << " ; j++) {"<<std::endl;
			  output << "dist_k" << i <<" += (data[global_index* " << dimensions
			  << "	+ j] - data[j + k_register" << i << " * " << dimensions << " ])"<<std::endl
			  <<"* (data[j + global_index* " << dimensions << " ] - data[j + k_register" << i << " * " << dimensions << " ]);"<<std::endl;
			  output << "}"<<std::endl;*/
		}
		return output.str();
	}
	std::string replace_max_k_register(size_t k) {
		std::stringstream output;
		output << this->indent[2] << this->floatType() << " tmp = dist;" <<std::endl;
		output << this->indent[2] << "int token = i;" <<std::endl;
		output << this->indent[2] << "int tmpi = i;" <<std::endl;
		for(size_t i = 0; i < k; i++) {
			output << this->indent[2] << "if(dist_k" << i <<" > dist) {"<<std::endl;
			output << this->indent[3] << "tmp = dist_k" << i <<";"<<std::endl;
			output << this->indent[3] << "dist_k" << i << " = dist;"<<std::endl;
			output << this->indent[3] << "dist = tmp;"<<std::endl;
			output << this->indent[3] << "tmpi = k_register" << i << ";"<<std::endl;
			output << this->indent[3] << "k_register" << i << " = token;"<<std::endl;
			output << this->indent[3] << "token = tmpi;"<<std::endl;
			output << this->indent[2] << "}"<<std::endl;
		}

		/*output << this->floatType() << " maxdist = 0.0;" <<std::endl;
		  for(size_t i = 0; i < k; i++) {
		  output << "if(dist_k" << i <<" > maxdist)"<<std::endl;
		  output << "maxdist = dist_k" << i << ";"<<std::endl;
		  }
		  output << "if(dist < maxdist) {"<<std::endl;
		  output << "if(dist_k0 == maxdist) {"<<std::endl;
		  output << "dist_k0 = dist;"<<std::endl;
		  output << "k_register0 = i;"<<std::endl;
		  output << "}"<<std::endl;
		  for(size_t i = 1; i < k; i++) {
		  output << "else if(dist_k" << i <<" == maxdist) {"<<std::endl;
		  output << "dist_k" << i <<" = dist;"<<std::endl;
		  output << "k_register" << i <<" = i;"<<std::endl;
		  output << "}"<<std::endl;
		  }*/
		//output << "}"<<std::endl;
		return output.str();
	}
	std::string copy_k_registers_into_global(size_t k) {
		std::stringstream output;
		for(size_t i = 0; i < k; i++) {
			output << this->indent[0] << "neighbors[global_index * "<< k <<" + " << i << "] = k_register" << i << ";" << std::endl;
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
				return this->reuseSource("DensityOCLMultiPlatform_create_graph.cl");
			}
		}

		std::stringstream sourceStream;

		if (this->floatType().compare("double") == 0) {
			sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
		}

		sourceStream<<""<<std::endl
					<<"void kernel connectNeighbors(global "<<this->floatType()<<" *data, global int *neighbors)" << std::endl
					<<"{"<<std::endl
					<< this->indent[0]	<< "size_t global_index = get_global_id(0);" << std::endl
					<< this->indent[0]	<< "size_t index = get_global_id(0);" << std::endl
					<< init_k_registers(k, dimensions)
					<< this->indent[0] << "for (unsigned int i = 0; i <	" << dataSize << "; i++) {" << std::endl
					<< this->indent[1] << "if (i != global_index) {"<<std::endl
					<<"//get distance to current point"<<std::endl
					<< this->indent[2] << this->floatType()<<" dist = 0.0;"<<std::endl
					<< this->indent[2] << "for (unsigned int j = 0; j <	 " << dimensions << " ; j++) {"<<std::endl
					<< this->indent[3] << "dist += (data[global_index* " << dimensions << "	 + j] - data[j + i* " << dimensions << " ])"<<std::endl
					<< this->indent[3] << "* (data[j + global_index* " << dimensions << " ] - data[j + i* " << dimensions << " ]);"<<std::endl
					<< this->indent[2] << "}" << std::endl
					<< replace_max_k_register(k)
					<< this->indent[1] << "}" << std::endl
					<< this->indent[0] << "}" << std::endl
					<< copy_k_registers_into_global(k)
					<< "}" << std::endl;
		if(kernelConfiguration.contains("WRITE_SOURCE")) {
			if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
				this->writeSource("DensityOCLMultiPlatform_create_graph.cl", sourceStream.str());
			}
		}


		return sourceStream.str();
	}

}
	;

}
}
}
