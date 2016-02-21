
#pragma once

#include <fstream>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename real_type>
class SourceBuilderMult: public base::KernelSourceBuilderBase<real_type> {
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

    SourceBuilderMult(std::shared_ptr<base::OCLDevice> device, json::Node &kernelConfiguration, size_t dims) :
            device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    }

    std::string generateSource(size_t dimensions, size_t gridsize) {
        std::stringstream sourceStream;

        if (this->floatType().compare("double") == 0) {
            sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl << std::endl;
			}

		sourceStream << "void kernel multdensity(global const int *starting_points,global const "<<this->floatType()<<" *alpha, global "<<this->floatType()<<" *result,const "<<this->floatType()<<" lambda)" << std::endl;
		sourceStream << "{" << std::endl;
		sourceStream << ""<<this->floatType()<<" gesamtint=0.0;" << std::endl;
		sourceStream << "for(private int i=0;i< " << gridsize << ";i++) {" << std::endl;
		sourceStream << ""<<this->floatType()<<" zellenintegral=1.0;" << std::endl;
		sourceStream << "for(private int dim=0;dim< " << dimensions << ";dim++) {" << std::endl;
		sourceStream << "int index=starting_points[get_global_id(0)* " << dimensions << "*2+2*dim];" << std::endl;
		sourceStream << "int level=starting_points[get_global_id(0)* " << dimensions << "*2+2*dim+1];" << std::endl;
		sourceStream << "int index2=starting_points[i* " << dimensions << "*2+2*dim];" << std::endl;
		sourceStream << "int level2=starting_points[i* " << dimensions << "*2+2*dim+1];" << std::endl;
		sourceStream << "if(starting_points[get_global_id(0)* " << dimensions << "*2+2*dim+1]>starting_points[i* " << dimensions << "*2+2*dim+1]) {" << std::endl;
		sourceStream << "index=starting_points[i* " << dimensions << "*2+2*dim];" << std::endl;
		sourceStream << "level=starting_points[i* " << dimensions << "*2+2*dim+1];" << std::endl;
		sourceStream << "index2=starting_points[get_global_id(0)* " << dimensions << "*2+2*dim];" << std::endl;
		sourceStream << "level2=starting_points[get_global_id(0)* " << dimensions << "*2+2*dim+1];" << std::endl;
		sourceStream << "}" << std::endl;
		sourceStream << "int teiler = (1 << level2);" << std::endl;
		sourceStream << ""<<this->floatType()<<" h=1.0 / teiler;" << std::endl;
		//sourceStream << "for(unsigned int z=1;z<=level2;z++)" << std::endl;
		//sourceStream << "h*=0.5;" << std::endl;
		sourceStream << ""<<this->floatType()<<" grenze1=h*(index2-1);" << std::endl;
		sourceStream << ""<<this->floatType()<<" grenze2=h*(index2+1);" << std::endl;
		sourceStream << "int u= (1 << level);" << std::endl;
		sourceStream << ""<<this->floatType()<<" uright=u*grenze2-index;" << std::endl;
		sourceStream << ""<<this->floatType()<<" uleft=u*grenze1-index;" << std::endl;
		sourceStream << "if(uleft<0)" << std::endl;
		sourceStream << "uleft*=-1;" << std::endl;
		sourceStream << "uleft=1-uleft;" << std::endl;
		sourceStream << "if(uleft<0)" << std::endl;
		sourceStream << "uleft=0;" << std::endl;
		sourceStream << "if(uright<0)" << std::endl;
		sourceStream << "uright*=-1;" << std::endl;
		sourceStream << "uright=1-uright;" << std::endl;
		sourceStream << "if(uright<0)" << std::endl;
		sourceStream << "uright=0;" << std::endl;
		sourceStream << "if(starting_points[i* " << dimensions << "*2+2*dim+1]==starting_points[get_global_id(0)* " << dimensions << "*2+2*dim+1]) {" << std::endl;
		sourceStream << "zellenintegral*=2.0/3.0*h;" << std::endl;
		sourceStream << "if(starting_points[i* " << dimensions << "*2+2*dim]!=starting_points[get_global_id(0)* " << dimensions << "*2+2*dim])" << std::endl;
		sourceStream << "zellenintegral=0.0;" << std::endl;
		sourceStream << "}" << std::endl;
		sourceStream << "else" << std::endl;
		sourceStream << "zellenintegral*=h/2.0*(uleft+uright);" << std::endl;
		sourceStream << "}" << std::endl;
		sourceStream << "gesamtint+=zellenintegral*alpha[i];" << std::endl;
		sourceStream << "}" << std::endl;
		sourceStream << "result[get_global_id(0)]=gesamtint;" << std::endl;
		sourceStream << "result[get_global_id(0)]+=alpha[get_global_id(0)]*lambda;" << std::endl;
		sourceStream << "}" << std::endl;

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
