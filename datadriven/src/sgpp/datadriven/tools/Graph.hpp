// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <sgpp/datadriven/tools/vpTree/VpTree.hpp>

#include <queue>

namespace sgpp {
namespace datadriven{

typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS,
boost::property<boost::vertex_index_t, size_t>, boost::no_property> UndirectedGraph;
/**
* @brief Class that encapsulates all methods and properties of a nearest neighbors graph
*/
class Graph {
 public:
  Graph(size_t vertices);

  Graph (const Graph &rhs) {
    this->graph = new UndirectedGraph(*(rhs.graph));
  }

  Graph &operator=(const Graph &rhs) {
    this->graph = new UndirectedGraph(*(rhs.graph));
    return *this;
  }

  ~Graph() = default;

  void addVertex();

  void disconnectVertex(size_t vertex);

  void createEdges(size_t vertex, std::priority_queue<VpHeapItem> nearestNeighbors);

  void addEdge(size_t vertex1, size_t vertex2);

  void deleteEdge(size_t vertex1, size_t vertex2);

  UndirectedGraph* getGraph();

  size_t getConnectedComponents(std::map<UndirectedGraph::vertex_descriptor, size_t> &componentMap);

 private:
    UndirectedGraph* graph;
};

}
}