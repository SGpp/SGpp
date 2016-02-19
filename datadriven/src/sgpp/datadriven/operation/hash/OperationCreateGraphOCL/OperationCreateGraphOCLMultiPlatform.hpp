
#pragma once

#include <omp.h>
#include <chrono>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/QueueLoadBalancer.hpp>
#include "OperationCreateGraphOCL.hpp"
#include "KernelCreateGraph.hpp"

namespace SGPP {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

template<typename T>
class OperationCreateGraphOCLMultiPlatform : public SGPP::datadriven::StreamingOCLMultiPlatform::OperationCreateGraphOCL
{
private:
	size_t dims;
	json::Node &configuration;
	SGPP::datadriven::DensityOCLMultiPlatform::KernelCreateGraph<T> *graph_kernel;
	std::vector<std::shared_ptr<base::OCLDevice>> devices;
	bool verbose;
	std::shared_ptr<base::OCLManagerMultiPlatform> manager;
	std::vector<T> dataVector;
public:

	OperationCreateGraphOCLMultiPlatform(base::DataMatrix& data, size_t dimensions,
									 std::shared_ptr<base::OCLManagerMultiPlatform> manager,
										 json::Node &configuration, size_t k) : OperationCreateGraphOCL(),dims(dimensions),configuration(configuration),devices(manager->getDevices()),manager(manager),
																				dataVector(data.getSize())
	{
		verbose=true;
		double *data_raw = data.getPointer();
		for(size_t i = 0; i < data.getSize(); i++)
			dataVector[i]=data_raw[i];
		graph_kernel=new SGPP::datadriven::DensityOCLMultiPlatform::KernelCreateGraph<T>(devices[0], dims, k, dataVector,
																						 manager, configuration);
	}

	~OperationCreateGraphOCLMultiPlatform() {
		delete graph_kernel;
	}

	void create_graph( std::vector<int> &resultVector)  {
		if(verbose)
			std::cout<<"Creating graph for"<<dataVector.size()<<" datapoints"<<std::endl;
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		try {
			this->graph_kernel->create_graph(resultVector);
		}
		catch(base::operation_exception &e) {
			std::cerr<<"Error! Could not create graph."<<std::endl
					 <<"Error Message: "<<e.what()<<std::endl;
			return;
		}
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;

		if (verbose)
		{
			std::cout << "duration create graph" << elapsed_seconds.count() << std::endl;
		}
	}
};

}
}
}
