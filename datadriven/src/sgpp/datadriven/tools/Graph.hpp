// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_BOOST_GRAPH

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <sgpp/datadriven/tools/vpTree/VpTree.hpp>

#include <map>
#include <queue>

namespace sgpp {
namespace datadriven {
/**
 * Definition of the type UndirecetdGraph
 */
typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS,
boost::property<boost::vertex_index_t, size_t>, boost::no_property> UndirectedGraph;

/**
* @brief Class that encapsulates all methods and properties of an unidrected Graph
*/
class Graph {
 public:
 /**
  * Constructs with no edges and a certan number of vertices
  * @param vertices Number of vertices contained in the graph
  */
  explicit Graph(size_t vertices);

  /**
   * Copy constructor
   * @param rhs copy object
   */
  Graph(const Graph &rhs) {
    this->graph = new UndirectedGraph(*(rhs.graph));
  }

  /**
   * Destructor
   */
  ~Graph() {
    delete this->graph;
  };

  /**
   * Method to add an additional vertex to the graph
   */
  void addVertex();

  /**
   * Removes a given vertex
   * @param vertex The index used to identify the vertex to remove
   */
  void removeVertex(size_t vertex);

  /**
   * Creates the edges of a vertex given its nearest neighbors
   * @param vertex The index used to identify the vertex to remove
   * @param nearestNeighbors priority queue obtaiend from a VP Tree with the indexes of the vertex
   * which are the nearest neighbors
   */
  void createEdges(size_t vertex, std::priority_queue<VpHeapItem> nearestNeighbors);

  /**
   * Adds an edge between to vertices
   * @param vertex1 The index used to identify the source vertex
   * @param vertex2 The index used to identify the sink vertex
   */
  void addEdge(size_t vertex1, size_t vertex2);

  /**
   * Deletes and edge between to vertices
   * @param vertex1
   * @param vertex2
   */
  void deleteEdge(size_t vertex1, size_t vertex2);

  /**
   * Gets the underlying boost graph structure
   * @return Pointer to the boost graph data structure
   */
  UndirectedGraph* getGraph();

  /**
   * Obtains the number of connected components and stores the labels in a given map
   * @param componentMap Map which contains the mapping from vertices to assigned labels
   * @return  Number of connected components detected
   */
  size_t getConnectedComponents(std::map<UndirectedGraph::vertex_descriptor, size_t> &componentMap);

 private:
    /**
     * Boost graph data structure
     */
    UndirectedGraph* graph;
};
}  // namespace datadriven
}  // namspace sgpp

#endif