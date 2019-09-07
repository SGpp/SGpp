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

#include <iostream>

namespace sgpp {
namespace datadriven {

class Cell {
  size_t dimension;
  double* corner;
  double* width;

 public:
  explicit Cell(size_t inp_dimension);
  Cell(size_t inp_dimension, double* inp_corner, double* inp_width);
  ~Cell();
  double getCorner(size_t d);
  double getWidth(size_t d);
  void setCorner(size_t d, double val);
  void setWidth(size_t d, double val);
  bool containsPoint(double point[]);
};


class SPTree {
  // Fixed constants
  static const size_t QT_NODE_CAPACITY = 1;

  // A buffer we use when doing force computations
  double* buff;

  // Properties of this node in the tree
  SPTree* parent;
  size_t dimension;
  bool is_leaf;
  size_t size;
  size_t cum_size;

  // Axis-aligned bounding box stored as a center with half-dimensions to
  // represent the boundaries of this quad tree
  Cell* boundary;

  // Indices in this space-partitioning tree node, corresponding center-of-mass,
  // and list of all children
  double* data;
  double* center_of_mass;
  size_t index[QT_NODE_CAPACITY];

  // Children
  SPTree** children;
  size_t no_children;

 public:
  SPTree(size_t D, double* inp_data, size_t N);
  SPTree(size_t D, double* inp_data, double* inp_corner, double* inp_width);
  SPTree(size_t D, double* inp_data, size_t N, double* inp_corner,
    double* inp_width);
  SPTree(SPTree* inp_parent, size_t D, double* inp_data, size_t N,
    double* inp_corner, double* inp_width);
  SPTree(SPTree* inp_parent, size_t D, double* inp_data,
    double* inp_corner, double* inp_width);
  ~SPTree();
  void setData(double* inp_data);
  SPTree* getParent();
  void construct(Cell boundary);
  bool insert(size_t new_index);
  void subdivide();
  bool isCorrect();
  void rebuildTree();
  void getAllIndices(size_t* indices);
  size_t getDepth();
  void computeNonEdgeForces(size_t point_index, double theta, double neg_f[], double* sum_Q);
  void computeEdgeForces(size_t* row_P, size_t* col_P,
    double* val_P, size_t N, double* pos_f);

 private:
  void init(SPTree* inp_parent, size_t D, double* inp_data, double*
    inp_corner, double* inp_width);
  void fill(size_t N);
  size_t getAllIndices(size_t* indices, size_t loc);
  bool isChild(size_t test_index, size_t start, size_t end);
};

}  // namespace datadriven
}  // namespace sgpp
