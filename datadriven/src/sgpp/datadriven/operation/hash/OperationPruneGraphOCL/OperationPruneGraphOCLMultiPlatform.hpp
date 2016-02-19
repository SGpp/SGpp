
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
#include "OperationPruneGraphOCL.hpp"
#include "KernelPruneGraph.hpp"

namespace SGPP {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

template<typename T>
class OperationPruneGraphOCLMultiPlatform : public SGPP::datadriven::StreamingOCLMultiPlatform::OperationPruneGraphOCL
{
private:
	size_t dims;
	size_t gridSize;
	size_t dataSize;
	json::Node &configuration;
	SGPP::datadriven::DensityOCLMultiPlatform::KernelPruneGraph<T> *graph_kernel;
	std::vector<std::shared_ptr<base::OCLDevice>> devices;
	bool verbose;
	std::shared_ptr<base::OCLManagerMultiPlatform> manager;

	std::vector<int> pointsVector;
	std::vector<T> alphaVector;
	std::vector<T> dataVector;
public:

	OperationPruneGraphOCLMultiPlatform( base::Grid& grid, base::DataVector& alpha, base::DataMatrix &data, size_t dims,
										 std::shared_ptr<base::OCLManagerMultiPlatform> manager,
										 json::Node &configuration, T treshold, size_t k) :
		OperationPruneGraphOCL(), dims(dims), gridSize(grid.getStorage()->size()), dataSize(data.getSize()), configuration(configuration), dataVector(data.getSize()),
		devices(manager->getDevices()), manager(manager)
	{
		verbose=true;
		//Store Grid in a opencl compatible buffer
		SGPP::base::GridStorage* gridStorage = grid.getStorage();
		size_t pointscount=0;
		for(int i=0;i<gridSize;i++)
		{
			SGPP::base::HashGridIndex *point=gridStorage->get(i);
			pointscount++;
			for(int d=0;d<dims;d++)
			{
				pointsVector.push_back(point->getIndex(d));
				pointsVector.push_back(point->getLevel(d));
			}
		}
		if(verbose)
			std::cout<<"Grid stored into integer array! Number of gridpoints: "<<pointscount<<std::endl;
		for(size_t i=0;i<gridSize;i++)
			alphaVector.push_back(alpha[i]);
		double *data_raw = data.getPointer();
		for(size_t i = 0; i < data.getSize(); i++)
			dataVector[i]=data_raw[i];
		graph_kernel=new SGPP::datadriven::DensityOCLMultiPlatform::KernelPruneGraph<T>(devices[0], dims, treshold, k, manager, configuration, pointsVector, alphaVector, dataVector);
	}

	~OperationPruneGraphOCLMultiPlatform() {
		delete graph_kernel;
	}

	virtual void prune_graph(std::vector<int> &graph) {
		if(verbose)
			std::cout<<"Pruning graph for"<<graph.size()<<" nodes"<<std::endl;
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		try {
			this->graph_kernel->prune_graph(pointsVector, alphaVector, dataVector, graph);
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
			std::cout << "duration prune graph: " << elapsed_seconds.count() << std::endl;
		}
	}
};

}
}
}
