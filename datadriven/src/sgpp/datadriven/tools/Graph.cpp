// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_BOOST_GRAPH
#include <sgpp/datadriven/tools/Graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <ctime>
#include <queue>
#include <map>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

namespace sgpp {
namespace datadriven {

Graph::Graph(size_t vertices) {
  this->graph = new UndirectedGraph(vertices);
}


void Graph::addVertex() {
  boost::add_vertex(*graph);
}

void Graph::removeVertex(size_t vertex) {
  auto vertex_descriptor = boost::vertex(vertex, *graph);
  boost::clear_vertex(vertex_descriptor, *graph);
  boost::remove_vertex(vertex_descriptor, *graph);
}

size_t Graph::getConnectedComponents(
    std::map<UndirectedGraph::vertex_descriptor, size_t> &componentMap) {

  /*
   * This index has to be created in order for the connected components to work
   * while using lists as a vertex conatiner
   */
  boost::property_map<UndirectedGraph, boost::vertex_index_t>::type
      index = get(boost::vertex_index, *graph);
  boost::graph_traits<UndirectedGraph>::vertex_iterator vi, vend;
  boost::graph_traits<UndirectedGraph>::vertices_size_type cnt = 0;
  for (boost::tie(vi, vend) = boost::vertices(*graph); vi != vend; ++vi) {
    boost::put(index, *vi, cnt++);
  }

  auto numberComponents = boost::connected_components(*graph,
      boost::make_assoc_property_map(componentMap));
  return numberComponents;
}

void Graph::createEdges(size_t vertex, std::priority_queue<VpHeapItem> nearestNeighbors) {
  while (!nearestNeighbors.empty()) {
    addEdge(vertex, nearestNeighbors.top().index);
    nearestNeighbors.pop();
  }
}

void Graph::addEdge(size_t vertex1, size_t vertex2) {
  boost::add_edge(boost::vertex(vertex1, *graph), boost::vertex(vertex2, *graph), *graph);
}

void Graph::deleteEdge(size_t vertex1, size_t vertex2) {
  boost::remove_edge(boost::vertex(vertex1, *graph), boost::vertex(vertex2, *graph), *graph);
}

UndirectedGraph* Graph::getGraph() {
  return graph;
}
}  // namespace datadriven
}  // namespace sgpp
#endif