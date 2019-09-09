/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */

/**
 * Code originally taken from https://lvdmaaten.github.io/tsne/
 * It has been modified in order to be adapted to the SG++ datamining
 * pipeline structure and has been parallelized
 */

#include "sptree.hpp"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>


namespace sgpp {
namespace datadriven {

// Constructs cell
Cell::Cell(size_t inp_dimension) {
  dimension = inp_dimension;
  corner = new double[dimension];
  width  = new double[dimension];
}

Cell::Cell(size_t inp_dimension, double* inp_corner, double* inp_width) {
  dimension = inp_dimension;
  corner = new double[dimension];
  width  = new double[dimension];
  for (size_t d = 0; d < dimension; d++) {
    setCorner(d, inp_corner[d]);
  }
  for (size_t d = 0; d < dimension; d++) {
    setWidth(d,  inp_width[d]);
  }
}

// Destructs cell
Cell::~Cell() {
  delete[] corner;
  delete[] width;
}

double Cell::getCorner(size_t d) {
  return corner[d];
}

double Cell::getWidth(size_t d) {
  return width[d];
}

void Cell::setCorner(size_t d, double val) {
  corner[d] = val;
}

void Cell::setWidth(size_t d, double val) {
  width[d] = val;
}

// Checks whether a point lies in a cell
bool Cell::containsPoint(double point[]) {
  for (size_t d = 0; d < dimension; d++) {
    if (corner[d] - width[d] > point[d]) {
      return false;
    }
    if (corner[d] + width[d] < point[d]) {
      return false;
    }
  }
  return true;
}


// Default constructor for SPTree -- build tree, too!
SPTree::SPTree(size_t D, double* inp_data, size_t N) {
  // Compute mean, width, and height of current map (boundaries of SPTree)
  size_t nD = 0;
  double* mean_Y = new double[D];
  double*  min_Y = new double[D];
  for (size_t d = 0; d < D; d++) {
    min_Y[d] =  DBL_MAX;
  }
  double*  max_Y = new double[D];
  for (size_t d = 0; d < D; d++)  {
    max_Y[d] = -DBL_MAX;
  }
  for (size_t n = 0; n < N; n++) {
    for (size_t d = 0; d < D; d++) {
      mean_Y[d] += inp_data[n * D + d];
      if (inp_data[nD + d] < min_Y[d]) {
        min_Y[d] = inp_data[nD + d];
      }
      if (inp_data[nD + d] > max_Y[d]) {
        max_Y[d] = inp_data[nD + d];
      }
    }
    nD += D;
  }

  for (size_t d = 0; d < D; d++) {
    mean_Y[d] /= static_cast<double> (N);
  }

  // Construct SPTree
  double* width = new double[D];

  for (size_t d = 0; d < D; d++) {
    width[d] = fmax(max_Y[d] - mean_Y[d], mean_Y[d] - min_Y[d]) + 1e-5;
  }
  init(NULL, D, inp_data, mean_Y, width);
  fill(N);

  // Clean up memory
  delete[] mean_Y;
  delete[] max_Y;
  delete[] min_Y;
  delete[] width;
}


// Constructor for SPTree with particular size and parent -- build the tree, too!
SPTree::SPTree(size_t D, double* inp_data, size_t N,
  double* inp_corner, double* inp_width) {
  init(NULL, D, inp_data, inp_corner, inp_width);
  fill(N);
}


// Constructor for SPTree with particular size (do not fill the tree)
SPTree::SPTree(size_t D, double* inp_data, double* inp_corner, double* inp_width) {
  init(NULL, D, inp_data, inp_corner, inp_width);
}


// Constructor for SPTree with particular size and parent (do not fill tree)
SPTree::SPTree(SPTree* inp_parent, size_t D, double* inp_data,
  double* inp_corner, double* inp_width) {
  init(inp_parent, D, inp_data, inp_corner, inp_width);
}


// Constructor for SPTree with particular size and parent -- build the tree, too!
SPTree::SPTree(SPTree* inp_parent, size_t D, double* inp_data,
  size_t N, double* inp_corner, double* inp_width) {
  init(inp_parent, D, inp_data, inp_corner, inp_width);
  fill(N);
}


// Main initialization function
void SPTree::init(SPTree* inp_parent, size_t D, double* inp_data,
  double* inp_corner, double* inp_width) {
  parent = inp_parent;
  dimension = D;
  no_children = 2;
  for (size_t d = 1; d < D; d++) {
    no_children *= 2;
  }
  data = inp_data;
  is_leaf = true;
  size = 0;
  cum_size = 0;

  boundary = new Cell(dimension);
  for (size_t d = 0; d < D; d++) {
    boundary->setCorner(d, inp_corner[d]);
  }
  for (size_t d = 0; d < D; d++) {
    boundary->setWidth(d, inp_width[d]);
  }
  children = new SPTree*[no_children];
  // children = reinterpret_cast<SPTree**> (malloc(no_children * sizeof(SPTree*)));
  for (size_t i = 0; i < no_children; i++) {
    children[i] = NULL;
  }

  center_of_mass = new double[D];;
  for (size_t d = 0; d < D; d++) {
    center_of_mass[d] = .0;
  }

  buff = new double[D];
}


// Destructor for SPTree
SPTree::~SPTree() {
  for (size_t i = 0; i < no_children; i++) {
    if (children[i] != NULL) delete children[i];
  }
  delete[] children;
  delete[] center_of_mass;
  delete[] buff;
  delete boundary;
}


// Update the data underlying this tree
void SPTree::setData(double* inp_data) {
  data = inp_data;
}


// Get the parent of the current tree
SPTree* SPTree::getParent() {
  return parent;
}


// Insert a point into the SPTree
bool SPTree::insert(size_t new_index) {
  // Ignore objects which do not belong in this quad tree
  double* point = data + new_index * dimension;
  if (!boundary->containsPoint(point)) {
      return false;
  }

  // Online update of cumulative size and center-of-mass
  cum_size++;
  double mult1 = static_cast<double> (cum_size - 1) / static_cast<double> (cum_size);
  double mult2 = 1.0 / static_cast<double> (cum_size);
  for (size_t d = 0; d < dimension; d++) {
    center_of_mass[d] *= mult1;
  }
  for (size_t d = 0; d < dimension; d++) {
    center_of_mass[d] += mult2 * point[d];
  }

  // If there is space in this quad tree and it is a leaf, add the object here
  if (is_leaf && size < QT_NODE_CAPACITY) {
    index[size] = new_index;
    size++;
    return true;
  }

  // Don't add duplicates for now (this is not very nice)
  bool any_duplicate = false;
  for (size_t n = 0; n < size; n++) {
    bool duplicate = true;
    for (size_t d = 0; d < dimension; d++) {
      if (point[d] != data[index[n] * dimension + d]) {
        duplicate = false;
        break;
      }
    }
    any_duplicate = any_duplicate | duplicate;
  }
  if (any_duplicate) {
    return true;
  }

  // Otherwise, we need to subdivide the current cell
  if (is_leaf) {
    subdivide();
  }

  // Find out where the point can be inserted
  for (size_t i = 0; i < no_children; i++) {
    if (children[i]->insert(new_index)) {
      return true;
    }
  }

  // Otherwise, the point cannot be inserted (this should never happen)
  return false;
}

// Create four children which fully divide this cell into four quads of equal area
void SPTree::subdivide() {
  // Create new children
  double* new_corner = new double[dimension];
  double* new_width  = new double[dimension];
  for (size_t i = 0; i < no_children; i++) {
    size_t div = 1;
    for (size_t d = 0; d < dimension; d++) {
      new_width[d] = .5 * boundary->getWidth(d);
      if ((i / div) % 2 == 1) {
        new_corner[d] = boundary->getCorner(d) - .5 * boundary->getWidth(d);
      } else {
        new_corner[d] = boundary->getCorner(d) + .5 * boundary->getWidth(d);
      }
      div *= 2;
    }
    children[i] = new SPTree(this, dimension, data, new_corner, new_width);
  }
  delete[] new_corner;
  delete[] new_width;

  // Move existing points to correct children
  for (size_t i = 0; i < size; i++) {
    bool success = false;
    for (size_t j = 0; j < no_children; j++) {
      if (!success) {
        success = children[j]->insert(index[i]);
      }
    }
    index[i] = -1;
  }

  // Empty parent node
  size = 0;
  is_leaf = false;
}

// Build SPTree on dataset
void SPTree::fill(size_t N) {
  for (size_t i = 0; i < N; i++) {
    insert(i);
  }
}


// Checks whether the specified tree is correct
bool SPTree::isCorrect() {
  for (size_t n = 0; n < size; n++) {
    double* point = data + index[n] * dimension;
    if (!boundary->containsPoint(point)) {
      return false;
    }
  }
  if (!is_leaf) {
    bool correct = true;
    for (size_t i = 0; i < no_children; i++) {
      correct = correct && children[i]->isCorrect();
    }
    return correct;
  } else {
    return true;
  }
}



// Build a list of all indices in SPTree
void SPTree::getAllIndices(size_t* indices) {
  getAllIndices(indices, 0);
}


// Build a list of all indices in SPTree
size_t SPTree::getAllIndices(size_t* indices, size_t loc) {
  // Gather indices in current quadrant
  for (size_t i = 0; i < size; i++) {
    indices[loc + i] = index[i];
  }
  loc += size;

  // Gather indices in children
  if (!is_leaf) {
    for (size_t i = 0; i < no_children; i++) {
      loc = children[i]->getAllIndices(indices, loc);
    }
  }
  return loc;
}


size_t SPTree::getDepth() {
  if (is_leaf) {
    return 1;
  }
  size_t depth = 0;
  for (size_t i = 0; i < no_children; i++) {
    depth = static_cast<size_t>(fmax(static_cast<double>(depth),
      static_cast<double>(children[i]->getDepth())));
  }
  return 1 + depth;
}


// Compute non-edge forces using Barnes-Hut algorithm
void SPTree::computeNonEdgeForces(size_t point_index, double theta,
  double neg_f[], double* sum_Q) {
  // Make sure that we spend no time on empty nodes or self-interactions
  if (cum_size == 0 || (is_leaf && size == 1 && index[0] == point_index)) {
    return;
  }

  // Compute distance between point and center-of-mass
  double D = .0;
  size_t ind = point_index * dimension;
  for (size_t d = 0; d < dimension; d++) {
    buff[d] = data[ind + d] - center_of_mass[d];
  }
  for (size_t d = 0; d < dimension; d++) {
    D += buff[d] * buff[d];
  }

  // Check whether we can use this node as a "summary"
  double max_width = 0.0;
  double cur_width;
  for (size_t d = 0; d < dimension; d++) {
      cur_width = boundary->getWidth(d);
      max_width = (max_width > cur_width) ? max_width : cur_width;
  }
  if (is_leaf || max_width / sqrt(D) < theta) {
      // Compute and add t-SNE force between point and current node
      D = 1.0 / (1.0 + D);
      double mult = static_cast<double>(cum_size) * static_cast<double>(D);
      *sum_Q += mult;
      mult *= D;
      for (size_t d = 0; d < dimension; d++) {
        neg_f[d] += mult * buff[d];
      }
  } else {
      // Recursively apply Barnes-Hut to children
      for (size_t i = 0; i < no_children; i++) {
        children[i]->computeNonEdgeForces(point_index, theta, neg_f, sum_Q);
      }
  }
}


// Computes edge forces
void SPTree::computeEdgeForces(size_t* row_P, size_t* col_P,
  double* val_P, size_t N, double* pos_f) {
  // Loop over all edges in the graph
  size_t ind1 = 0;
  size_t ind2 = 0;
  double D;
  for (size_t n = 0; n < N; n++) {
    for (size_t i = row_P[n]; i < row_P[n + 1]; i++) {
      // Compute pairwise distance and Q-value
      D = 1.0;
      ind2 = col_P[i] * dimension;
      for (size_t d = 0; d < dimension; d++) {
        buff[d] = data[ind1 + d] - data[ind2 + d];
      }
      for (size_t d = 0; d < dimension; d++) {
        D += buff[d] * buff[d];
      }
      D = val_P[i] / D;

      // Sum positive force
      for (size_t d = 0; d < dimension; d++) {
        pos_f[ind1 + d] += D * buff[d];
      }
    }
    ind1 += dimension;
}
}

}  // namespace datadriven
}  // namespace sgpp
