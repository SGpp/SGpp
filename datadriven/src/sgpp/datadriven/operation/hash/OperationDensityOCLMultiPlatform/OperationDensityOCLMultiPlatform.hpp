
// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
#include "OperationDensityOCL.hpp"
#include "KernelMult.hpp"
#include "KernelB.hpp"

namespace SGPP {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

template<typename T>
class OperationDensityOCLMultiPlatform: public SGPP::datadriven::StreamingOCLMultiPlatform::OperationDensityOCL
{
private:
	size_t dims;
	size_t gridSize;
	SGPP::base::Grid& grid;
	json::Node &firstKernelConfig;
	json::Node &secondKernelConfig;
	SGPP::datadriven::DensityOCLMultiPlatform::KernelDensityMult<T> *multKernel;
	SGPP::datadriven::DensityOCLMultiPlatform::KernelDensityB<T> *bKernel;
	std::vector<std::shared_ptr<base::OCLDevice>> devices;
	bool verbose;
	std::vector<int> points;
	std::shared_ptr<base::OCLManagerMultiPlatform> manager;
	T lambda;
public:

	OperationDensityOCLMultiPlatform(base::Grid& grid,size_t dimensions,
									 std::shared_ptr<base::OCLManagerMultiPlatform> manager,
									 json::Node &firstKernelConfig, json::Node &secondKernelConfig, T lambda) : OperationDensityOCL(),dims(dimensions),gridSize(grid.getStorage()->size()),grid(grid),firstKernelConfig(firstKernelConfig),
																												secondKernelConfig(secondKernelConfig), devices(manager->getDevices()),manager(manager),lambda(lambda)
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
				points.push_back(point->getIndex(d));
				points.push_back(point->getLevel(d));
			}
		}
		if(verbose)
			std::cout<<"Grid stored into integer array! Number of gridpoints: "<<pointscount<<std::endl;
		multKernel=new SGPP::datadriven::DensityOCLMultiPlatform::KernelDensityMult<T>(devices[0], dims,
																					   manager, firstKernelConfig, points, lambda);
		bKernel=new SGPP::datadriven::DensityOCLMultiPlatform::KernelDensityB<T>(devices[0], dims, manager, secondKernelConfig,
																				 points);
	}

	~OperationDensityOCLMultiPlatform() {
		delete multKernel;
		delete bKernel;
	}

	void mult(base::DataVector& alpha, base::DataVector& result) override {
		std::vector<T> alphaVector(gridSize);
		std::vector<T> resultVector(gridSize);
		for(size_t i=0;i<gridSize;i++)
		{
			alphaVector[i]=alpha[i];
			resultVector[i]=result[i];
		}
		if(verbose)
			std::cout<<"starting multiplication with "<<gridSize<<" entries"<<std::endl;
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
			this->multKernel->mult(alphaVector, resultVector);
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;

		if (verbose)
		{
			std::cout << "duration mult ocl: " << elapsed_seconds.count() << std::endl;
		}
		for(size_t i=0; i< gridSize; i++)
			result[i] = resultVector[i];
	}

	void generateb(base::DataMatrix &dataset, SGPP::base::DataVector &b) {
		std::chrono::time_point<std::chrono::system_clock> start, end;
		if (verbose) {
			std::cout << "starting rhs kernel methode! datasize: " <<b.getSize()<< std::endl;
		}
		std::vector<T> bVector(b.getSize());
		std::vector<T> datasetVector(dataset.getSize());
		double *data_raw = dataset.getPointer();
		for(size_t i = 0; i < dataset.getSize(); i++)
			datasetVector[i]=data_raw[i];
		start = std::chrono::system_clock::now();
		bKernel->rhs(datasetVector, bVector);
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;

		for(size_t i = 0; i < b.getSize(); i++)
			b[i]=bVector[i];
		if (verbose) {
			std::cout << "duration rhs ocl: " << elapsed_seconds.count() << std::endl;
		}
	}
};

}
}
}
