

// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>

namespace SGPP {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

class OperationCreateGraphOCL
{
protected:
	static int find_neighbors(int index, std::vector<int> &nodes, int cluster, int k, std::vector<int> &clusterList,
							  bool overwrite=false) {
		int currIndex;
		clusterList[index]=cluster;
		bool overwrite_enabled = overwrite;
		for (int i = index*k; i < (index+1) *k; i++) {
			if(nodes[i]==-2)
				continue;
			currIndex = nodes[i];
			if (nodes[currIndex*k]!=-1) {
				if(clusterList[currIndex]==0 && !overwrite_enabled) {
					clusterList[currIndex] = cluster;
					int ret_cluster = OperationCreateGraphOCL::find_neighbors(currIndex, nodes, cluster, k, clusterList);
					if(ret_cluster != cluster) {
						cluster = ret_cluster;
						clusterList[index]=cluster;
						overwrite_enabled = true;
						i = index*k;
						continue;
					}
				}
				else if(!overwrite_enabled && clusterList[currIndex]!=cluster) {
					cluster = clusterList[currIndex];
					clusterList[index]=cluster;
					overwrite_enabled = true;
					i = index*k-1;
					continue;
				}
				else{
					if(clusterList[currIndex]!=cluster) {
						OperationCreateGraphOCL::find_neighbors(currIndex, nodes, cluster, k, clusterList,true);
						clusterList[currIndex] = cluster;
					}
				}
			}
		}
		return cluster;
	}
public:
	OperationCreateGraphOCL()  {
	}

	virtual void create_graph(base::DataVector& data, std::vector<int> &resultVector) {}
	static std::vector<int> find_clusters(std::vector<int> &graph, size_t k) {
		std::vector<int> clusters(graph.size()/k);
		int clustercount = 0;
		for (size_t node = 0; node < clusters.size(); node++)
			clusters[node]=0;
		for (int i = 0; i < clusters.size(); i++) {
			if(clusters[i] == 0 && graph[i*k]!=-1) {
				clustercount++;
				if(OperationCreateGraphOCL::find_neighbors(i, graph, clustercount, k, clusters) != clustercount)
					clustercount--;
			}
		}
		std::cout<<"Found "<<clustercount<<" clusters!"<<std::endl;
		return clusters;

	}

};

}
}
}
