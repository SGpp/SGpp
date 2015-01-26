// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <list>
#include <vector>
using namespace std;
#include <iostream>

class Graph
{
    int V;    // No. of vertices
    vector<list<int> > *adj = NULL;    // An array of adjacency lists
    void dfs(int v, std::vector<int> &components, int label);
    std::vector<int> *components = NULL;
public:
    Graph(int V);
    ~Graph();
    void addEdge(int v, int w);
    std::vector<int> getComponents();
};

Graph::Graph(int V)
{
    this->V = V;
    adj = new vector<list<int> >(V);
}

Graph::~Graph()
{
	delete components;
	delete adj;
}

void Graph::addEdge(int v, int w)
{

    adj->at(v).push_back(w); // Add w to v’s list.
    adj->at(v).push_back(v);
}

void Graph::dfs(int v, std::vector<int> &components, int label)
{
    // Mark the current node as visited and print it
    components[v] = label;
    // Recur for all the vertices adjacent to this vertex
    list<int>::iterator i;
    for(i = (adj->at(v)).begin(); i != (adj->at(v)).end(); ++i){
        if(components[*i] == -1)
            dfs(*i, components, label);
    }
}


std::vector<int> Graph::getComponents()
{
    components = new std::vector<int>(V);

    // Mark all the vertices as not visited
    for(int i = 0; i < V; i++)
      (*components)[i] = -1;

    // Fill vertices in stack according to their finishing times
    int label = 0;
    for(int i = 0; i < V; i++){
        if((*components)[i] == -1){
            dfs(i, *components, label);
            label++;
        }
    }
    return *components;
}

