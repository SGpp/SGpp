// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#pragma once
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/KernelSourceBuilderBase.hpp>

#include <fstream>
#include <string>

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
            output << "ptrData[(" << dataBlockSize << " * globalIdx) + (resultSize * " << dim
                   << ") + " << dataBlockingIndex << "]";
        } else {
            std::string error("OCL Error: Illegal value for parameter \"KERNEL_STORE_DATA\"");
            throw new base::operation_exception(error.c_str());
        }
        return output.str();
    }

 public:
    SourceBuilderPruneGraph(std::shared_ptr<base::OCLDevice> device,
                            json::Node &kernelConfiguration, size_t dims) :
        device(device), kernelConfiguration(kernelConfiguration), dims(dims) {
    }

    std::string generateSource(size_t dimensions, size_t gridSize, size_t k, real_type treshold) {
        if (kernelConfiguration.contains("REUSE_SOURCE")) {
            if (kernelConfiguration["REUSE_SOURCE"].getBool()) {
                return this->reuseSource("DensityOCLMultiPlatform_prune_graph.cl");
            }
        }

        std::stringstream sourceStream;

        if (this->floatType().compare("double") == 0) {
            sourceStream << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable"
                         << std::endl << std::endl;
        }

        sourceStream << "" << std::endl
                     << "" << this->floatType() << " get_u(private const "
                     << this->floatType() << " grenze,private const int index," << std::endl
                     << "private const int level)" << std::endl
                     << "{" << std::endl
                     << this->indent[0] << "private " << this->floatType()
                     << " ret = (1 << level);" << std::endl
                     << this->indent[0] << "ret*=grenze;" << std::endl
                     << this->indent[0] << "ret-=index;" << std::endl
                     << this->indent[0] << "if (ret<0.0)" << std::endl
                     << this->indent[0] << "ret*=-1.0;" << std::endl
                     << this->indent[0] << "ret=1-ret;" << std::endl
                     << this->indent[0] << "if (ret<0.0)" << std::endl
                     << this->indent[1] << "ret=0.0;" << std::endl
                     << this->indent[0] << "return ret;" << std::endl
                     << "}" << std::endl
                     << "" << std::endl
                     << "void kernel removeEdges(global int *nodes,"
                     << "global const int *starting_points,global const " << this->floatType()
                     << " *data," << std::endl
                     << this->indent[0] << "global const " << this->floatType() << " *alphas) {"
                     << std::endl
                     << this->indent[0] << "size_t index = get_global_id(0);" << std::endl
                     << this->indent[0] << "size_t global_index=get_global_id(0);" << std::endl
                     << this->indent[0] << "" << this->floatType() << " endwert=0;" << std::endl
                     << this->indent[0] << "for (int gridpoint=0;gridpoint< " << gridSize
                     << " ;gridpoint++)" << std::endl
                     << this->indent[0] << "{" << std::endl
                     << this->indent[1] << "" << this->floatType() << " wert=1;" << std::endl
                     << this->indent[1] << "for (int dimension=0;dimension< " << dimensions
                     << ";dimension++)" << std::endl
                     << this->indent[1] << "{" << std::endl
                     << this->indent[2] << "wert*=get_u(data[global_index* " << dimensions
                     << "+dimension],starting_points[gridpoint*2* " << dimensions
                     << "+2*dimension]," << std::endl
                     << this->indent[3] << "starting_points[gridpoint*2* " << dimensions
                     << "+2*dimension+1]);" << std::endl
                     << this->indent[1] << "}" << std::endl
                     << this->indent[1] << "endwert+=wert*alphas[gridpoint];" << std::endl
                     << this->indent[0] << "}" << std::endl
                     << this->indent[0] << "if (endwert<0.0)" << std::endl
                     << this->indent[1] << "endwert=0.0;" << std::endl
                     << this->indent[0] << "if (endwert< " << treshold << " )" << std::endl
                     << this->indent[0] << "{" << std::endl
                     << this->indent[1] << "for (int i = 0; i <  " << k << " ; i++)" << std::endl
                     << this->indent[1] << "{" << std::endl
                     << this->indent[2] << "nodes[ " << k << " *index + i] = -1;" << std::endl
                     << this->indent[1] << "}" << std::endl
                     << this->indent[0] << "}" << std::endl
                     << this->indent[0] << "else //Remove Edges" << std::endl
                     << this->indent[0] << "{" << std::endl
                     << this->indent[1] << "for (int i = 0; i <  " << k << " ; i++)" << std::endl
                     << this->indent[1] << "{" << std::endl
                     << this->indent[2] << "//Calculate density" << std::endl
                     << this->indent[2] << "" << this->floatType() << " endwert=0;" << std::endl
                     << this->indent[2] << "int nachbar=nodes[index* " << k << " +i];" << std::endl
                     << this->indent[2] << "for (int gridpoint=0;gridpoint< " << gridSize
                     << " ;gridpoint++)" << std::endl
                     << this->indent[2] << "{" << std::endl
                     << this->indent[3] << "" << this->floatType() << " wert=1;" << std::endl
                     << this->indent[3] << "for (int dimension=0;dimension< " << dimensions
                     << ";dimension++)" << std::endl
                     << this->indent[3] << "{" << std::endl
                     << this->indent[4] << "" << this->floatType() << " dimension_point=0;"
                     << std::endl
                     << this->indent[4] << "if (data[global_index* " << dimensions
                     << "+dimension]>data[dimension+nachbar* " << dimensions << "])" << std::endl
                     << this->indent[5] << "dimension_point=data[dimension+nachbar* "
                     << dimensions << "]+" << std::endl
                     << this->indent[5] << "(data[global_index* " << dimensions
                     << "+dimension]-data[dimension+nachbar* " << dimensions << "])*0.5;"
                     << std::endl
                     << this->indent[4] << "else" << std::endl
                     << this->indent[5] << "dimension_point=data[global_index* " << dimensions
                     << "+dimension]+" << std::endl
                     << this->indent[5] << "(data[dimension+nachbar* " << dimensions
                     << "]-data[global_index* " << dimensions << "+dimension])*0.5;" << std::endl
                     << this->indent[4] << "wert*=get_u(dimension_point,"
                     << "starting_points[gridpoint*2* " << dimensions << "+2*dimension],"
                     << std::endl
                     << this->indent[4] << "starting_points[gridpoint*2* " << dimensions
                     << "+2*dimension+1]);" << std::endl
                     << this->indent[3] << "}" << std::endl
                     << this->indent[3] << "endwert+=wert*alphas[gridpoint];" << std::endl
                     << this->indent[2] << "}" << std::endl
                     << this->indent[2] << "if (endwert<0.0)" << std::endl
                     << this->indent[3] << "endwert=0.0;" << std::endl
                     << this->indent[2] << "if (endwert< " << treshold << " )" << std::endl
                     << this->indent[2] << "{" << std::endl
                     << this->indent[3] << "nodes[ " << k << " *index + i] = -2;" << std::endl
                     << this->indent[2] << "}" << std::endl
                     << this->indent[1] << "}" << std::endl
                     << this->indent[0] << "}" << std::endl
                     << "}" << std::endl;
        if (kernelConfiguration.contains("WRITE_SOURCE")) {
            if (kernelConfiguration["WRITE_SOURCE"].getBool()) {
                this->writeSource("DensityOCLMultiPlatform_prune_graph.cl", sourceStream.str());
            }
        }
        return sourceStream.str();
    }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace SGPP
