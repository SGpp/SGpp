// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONCREATEGRAPHOCL_H
#define OPERATIONCREATEGRAPHOCL_H

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <vector>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/KernelCreateGraph.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// Pure virtual base class for the k nearest neighbor opencl operation
class OperationCreateGraphOCL {
 protected:
  /// Recursive function for traversing the k nearest neighbor graph
  static size_t find_neighbors(size_t index, std::vector<int> &nodes, size_t cluster, size_t k,
                               std::vector<size_t> &clusterList,
                               bool overwrite = false) {
    size_t currIndex;
    clusterList[index] = cluster;
    bool overwrite_enabled = overwrite;
    bool removed = true;
    for (size_t i = index*k; i < (index+1) *k; i++) {
      if (nodes[i] == -2)
        continue;
      removed = false;
      currIndex = nodes[i];
      if (nodes[currIndex*k] != -1) {
        if (clusterList[currIndex] == 0 && !overwrite_enabled) {
          clusterList[currIndex] = cluster;
          size_t ret_cluster = OperationCreateGraphOCL::find_neighbors(currIndex, nodes,
                                                                       cluster, k,
                                                                       clusterList);
          if (ret_cluster != cluster) {
            cluster = ret_cluster;
            clusterList[index] = cluster;
            overwrite_enabled = true;
            i = index * k;
            continue;
          }
        } else if (!overwrite_enabled && clusterList[currIndex] != cluster) {
          cluster = clusterList[currIndex];
          clusterList[index] = cluster;
          overwrite_enabled = true;
          i = index * k - 1;
          continue;
        } else {
          if (clusterList[currIndex] != cluster) {
            OperationCreateGraphOCL::find_neighbors(currIndex, nodes, cluster, k,
                                                    clusterList, true);
            clusterList[currIndex] = cluster;
          }
        }
      }
    }
    if (removed) {
      clusterList[index] = 0;
      cluster = 0;
    }
    return cluster;
  }

 public:
  OperationCreateGraphOCL()  {
  }

  /// Pure virtual function to create the k nearest neighbor graph for some datapoints of a dataset
  virtual void create_graph(std::vector<int> &resultVector, int startid = 0,
                            int chunksize = 0) = 0;
  virtual void begin_graph_creation(int startid, int chunksize) = 0;
  virtual void finalize_graph_creation(std::vector<int> &resultVector, int startid,
                                       int chunksize) = 0;
  /// Assign a clusterindex for each datapoint using the connected components of the graph
  static std::vector<size_t> find_clusters(std::vector<int> &graph, size_t k) {
    std::vector<size_t> clusters(graph.size()/k);
    size_t clustercount = 0;
    for (size_t node = 0; node < clusters.size(); node++)
      clusters[node] = 0;
    for (size_t i = 0; i < clusters.size(); i++) {
      if (clusters[i] == 0 && graph[i*k] != -1) {
        clustercount++;
        if (OperationCreateGraphOCL::find_neighbors(i, graph, clustercount, k, clusters) !=
            clustercount)
          clustercount--;
      }
    }
    std::cout << "Found " << clustercount << " clusters!" << std::endl;
    return clusters;
  }

  virtual ~OperationCreateGraphOCL(void) {}
  /// Add the default parameters to the the configuration
  static void load_default_parameters(base::OCLOperationConfiguration *parameters) {
    if (parameters->contains("INTERNAL_PRECISION") == false) {
      std::cout << "Warning! No internal precision setting detected."
                << " Using double precision from now on!" << std::endl;
      parameters->addIDAttr("INTERNAL_PRECISION", "double");
    }
    if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
      DensityOCLMultiPlatform::KernelCreateGraph<float>::augmentDefaultParameters(*parameters);
    } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
      DensityOCLMultiPlatform::KernelCreateGraph<double>::augmentDefaultParameters(*parameters);
    } else {
      std::stringstream errorString;
      errorString << "Error creating operation\"CreateGraphOCL\": "
                  << " invalid value for parameter \"INTERNAL_PRECISION\"";
      throw base::operation_exception(errorString.str().c_str());
    }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp


#endif /* OPERATIONCREATEGRAPHOCL_H */
