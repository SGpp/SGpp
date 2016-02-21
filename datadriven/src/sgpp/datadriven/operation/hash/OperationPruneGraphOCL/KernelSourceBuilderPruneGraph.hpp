
#pragma once

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename real_type>
class SourceBuilderPruneGraph: public base::KernelSourceBuilderBase<real_type> {
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

	SourceBuilderPruneGraph(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration, size_t dims) :
		device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
	}

	std::string generateSource() {
		if(kernelConfiguration.contains("REUSE_SOURCE")) {
			if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
				return this->reuseSource("DensityOCLMultiPlatform_prune_graph.cl");
			}
		}

		std::stringstream sourceStream;

		if (this->floatType().compare("double") == 0) {
			sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
		}

		sourceStream<<""<<std::endl
					<<""<<this->floatType()<<" get_u(private const "<<this->floatType()<<" grenze,private const int index,"<<std::endl
					<<"		 private const int level)"<<std::endl
					<<"{"<<std::endl
					<<"private "<<this->floatType()<<" ret=1.0;"<<std::endl
					<<"for(int z=1;z<=level;z++)"<<std::endl
					<<"	ret*=2.0;"<<std::endl
					<<"ret*=grenze;"<<std::endl
					<<"ret-=index;"<<std::endl
					<<"if(ret<0.0)"<<std::endl
					<<"	ret*=-1.0;"<<std::endl
					<<"ret=1-ret;"<<std::endl
					<<"if(ret<0.0)"<<std::endl
					<<"	ret=0.0;"<<std::endl
					<<"return ret;"<<std::endl
					<<"}"<<std::endl
					<<""<<std::endl
					<<" void kernel removeEdges(global int *nodes,global const int *starting_points,global const "<<this->floatType()<<" *data,"<<std::endl
					<<"					global const "<<this->floatType()<<" *alphas,const int dimensions,const int gridsize,"<<std::endl
					<<"				const int k,const "<<this->floatType()<<" treshold) {"<<std::endl
					<<"size_t index = get_global_id(0);"<<std::endl
					<<"size_t global_index=get_global_id(0);"<<std::endl
					<<""<<this->floatType()<<" endwert=0;"<<std::endl
					<<"for(int gridpoint=0;gridpoint<gridsize;gridpoint++)"<<std::endl
					<<"{"<<std::endl
					<<"	"<<this->floatType()<<" wert=1;"<<std::endl
					<<"	for(int dimension=0;dimension<dimensions;dimension++)"<<std::endl
					<<"	{"<<std::endl
					<<"		wert*=get_u(data[global_index*dimensions+dimension],starting_points[gridpoint*2*dimensions+2*dimension],"<<std::endl
					<<"					starting_points[gridpoint*2*dimensions+2*dimension+1]);"<<std::endl
					<<"	}"<<std::endl
					<<"	endwert+=wert*alphas[gridpoint];"<<std::endl
					<<"}"<<std::endl
					<<"if(endwert<0.0)"<<std::endl
					<<"	endwert=0.0;"<<std::endl
					<<"if (endwert<treshold)"<<std::endl
					<<"{"<<std::endl
					<<"	for (int i = 0; i < k; i++)"<<std::endl
					<<"	{"<<std::endl
					<<"		nodes[k*index + i] = -1;"<<std::endl
					<<"	}"<<std::endl
					<<"}"<<std::endl
					<<"else //Remove Edges"<<std::endl
					<<"{"<<std::endl
					<<"	for (int i = 0; i < k; i++)"<<std::endl
					<<"	{"<<std::endl
					<<"		//Calculate density"<<std::endl
					<<"		"<<this->floatType()<<" endwert=0;"<<std::endl
					<<"	int nachbar=nodes[index*k+i];"<<std::endl
					<<"		for(int gridpoint=0;gridpoint<gridsize;gridpoint++)"<<std::endl
					<<"		{"<<std::endl
					<<"			"<<this->floatType()<<" wert=1;"<<std::endl
					<<"			for(int dimension=0;dimension<dimensions;dimension++)"<<std::endl
					<<"			{"<<std::endl
					<<"				"<<this->floatType()<<" dimension_point=0;"<<std::endl
					<<"				if(data[global_index*dimensions+dimension]>data[dimension+nachbar*dimensions])"<<std::endl
					<<"					dimension_point=data[dimension+nachbar*dimensions]+"<<std::endl
					<<"						(data[global_index*dimensions+dimension]-data[dimension+nachbar*dimensions])*0.5;"<<std::endl
					<<"				else"<<std::endl
					<<"					dimension_point=data[global_index*dimensions+dimension]+"<<std::endl
					<<"						(data[dimension+nachbar*dimensions]-data[global_index*dimensions+dimension])*0.5;"<<std::endl
					<<"				wert*=get_u(dimension_point,starting_points[gridpoint*2*dimensions+2*dimension],"<<std::endl
					<<"							starting_points[gridpoint*2*dimensions+2*dimension+1]);"<<std::endl
					<<"			}"<<std::endl
					<<"			endwert+=wert*alphas[gridpoint];"<<std::endl
					<<"		}"<<std::endl
					<<"		if(endwert<0.0)"<<std::endl
					<<"			endwert=0.0;"<<std::endl
					<<"		if (endwert<treshold)"<<std::endl
					<<"		{"<<std::endl
					<<"			nodes[k*index + i] = -2;"<<std::endl
					<<"		}"<<std::endl
					<<"	}"<<std::endl
					<<"}"<<std::endl
					<<"}"<<std::endl;
		if(kernelConfiguration.contains("WRITE_SOURCE")) {
			if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
				this->writeSource("DensityOCLMultiPlatform_prune_graph.cl", sourceStream.str());
			}
		}


		return sourceStream.str();
	}

}
	;

}
}
}
