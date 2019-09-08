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
 *
 */

/**
 * Code originally taken from https://lvdmaaten.github.io/tsne/
 * It has been modified in order to be adapted to the SG++ datamining
 * pipeline structure and has been parallelized
 */

#include <sgpp/datadriven/datamining/modules/visualization/algorithms/bhtsne/tsne.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/algorithms/bhtsne/vptree.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/algorithms/bhtsne/sptree.hpp>
#include <omp.h>
#include <cfloat>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>

namespace sgpp {
namespace datadriven {

TSNE::TSNE() {
}

// Perform t-SNE
void TSNE::run(double* X, size_t N, size_t D, double* Y, size_t no_dims, double perplexity,
  double theta, size_t rand_seed, bool skip_random_init, size_t max_iter, size_t mom_switch_iter) {

  srand(static_cast<unsigned int>(rand_seed));
  // Determine whether we are using an exact algorithm
  if (static_cast<double>(N - 1) < 3 * perplexity) {
    printf("Perplexity too large for the number of data points!\n");
    exit(1);
  }
  printf("Using no_dims = %zu, perplexity = %f, and theta = %f\n", no_dims, perplexity, theta);
  bool exact = (theta == .0) ? true : false;

  // Set learning parameters
  float total_time = .0;
  clock_t start, end;
  double momentum = .5, final_momentum = .8;
  double eta = 200.0;

  // Allocate some memory
  static std::unique_ptr<double[]> dY (new double[N * no_dims]);
  static std::unique_ptr<double[]> uY  (new double[N * no_dims]);
  static std::unique_ptr<double[]> gains (new double[N * no_dims]);

  for (size_t i = 0; i < N * no_dims; i++) {
    uY[i] =  .0;
  }
  for (size_t i = 0; i < N * no_dims; i++) {
    gains[i] = 1.0;
  }

  // Normalize input data (to prevent numerical problems)
  printf("Computing input similarities...\n");
  start = clock();
  zeroMean(X, N, D);
  double max_X = .0;

  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < N * D; i++) {
    if (fabs(X[i]) > max_X) {
      max_X = fabs(X[i]);
    }
  }
  for (size_t i = 0; i < N * D; i++) {
    X[i] /= max_X;
  }

  // Compute input similarities for exact t-SNE

  // Allocate the memory we need
  size_t K = static_cast<size_t> (3 * perplexity);
  static std::unique_ptr<size_t[]> row_P (new size_t[N + 1]);
  static std::unique_ptr<size_t[]> col_P (new size_t[N*K]);
  static std::unique_ptr<double[]> val_P (new double[N*K]);

  // Compute input similarities for approximate t-SNE
  // Compute asymmetric pairwise input similarities
  computeGaussianPerplexity(X, N, D, row_P.get(), col_P.get(), val_P.get(), perplexity,
    K);

  // Symmetrize input similarities
  symmetrizeMatrix(row_P, col_P, val_P, N);
  double sum_P = .0;
  for (size_t i = 0; i < row_P[N]; i++) {
    sum_P += val_P[i];
  }

  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < row_P[N]; i++) {
    val_P[i] /= sum_P;
  }

  end = clock();

  // Initialize solution (randomly)
  if (skip_random_init != true) {
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < N * no_dims; i++) {
      Y[i] = randn() * .0001;
    }
  }

  printf("Input similarities computed in %4.2f seconds "
    "(sparsity = %f)!\nLearning embedding...\n",
    static_cast<double> ((end - start) / CLOCKS_PER_SEC),
    static_cast<double> (row_P[N]) /
    (static_cast<double> (N) * static_cast<double> (N)));

  start = clock();
  for (size_t iter = 0; iter < max_iter; iter++) {

    computeGradient(row_P.get(), col_P.get(), val_P.get(), Y, N, no_dims, dY.get(), theta);
    // Update gains
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < N * no_dims; i++) {
      gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
      if (gains[i] < .01) {
        gains[i] = .01;
      }
      // Perform gradient update (with momentum and gains)
      uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
      Y[i] = Y[i] + uY[i];
    }
    // Make solution zero-mean
    zeroMean(Y, N, no_dims);

    if (iter == mom_switch_iter) {
      momentum = final_momentum;
    }

    // Print out progress
    if (iter > 0 && (iter % 50 == 0 || iter == max_iter - 1)) {
      end = clock();
      double C = .0;
        // doing approximate computation here!
      C = evaluateError(row_P.get(), col_P.get(), val_P.get(), Y, N, no_dims, theta);
      if (iter == 0) {
          printf("Iteration %zu: error is %f\n", iter + 1, C);
      } else {
        total_time += static_cast<float> ((end - start) / CLOCKS_PER_SEC);
        printf("Iteration %zu: error is %f (50 iterations in %4.2f seconds)\n",
          iter, C, static_cast<float> ((end - start) / CLOCKS_PER_SEC));
      }
      start = clock();
    }
  }
  end = clock(); total_time += static_cast<float> ((end - start) / CLOCKS_PER_SEC);

  // Clean up memory
  printf("Fitting performed in %4.2f seconds.\n", total_time);
}


// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
void TSNE::computeGradient(size_t* inp_row_P,
  size_t* inp_col_P, double* inp_val_P, double* Y,
  size_t N, size_t D, double* dC, double theta) {
  // Construct space-partitioning tree on current map
  SPTree* tree = new SPTree(D, Y, N);

  // Compute all terms required for t-SNE gradient
  double sum_Q = .0;
  double* pos_f = reinterpret_cast<double*> (calloc(N * D, sizeof(double)));
  double* neg_f = reinterpret_cast<double*> (calloc(N * D, sizeof(double)));
  if (pos_f == NULL || neg_f == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }

  #pragma omp parallel sections
  {
    #pragma omp section
    {
      tree->computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);
    }

    #pragma omp section
    {
      #pragma omp parallel for schedule(dynamic)
      for (size_t n = 0; n < N; n++) {
        tree->computeNonEdgeForces(n, theta, neg_f + n * D, &sum_Q);
      }
    }
  }

  // Compute final t-SNE gradient
  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < N * D; i++) {
    dC[i] = pos_f[i] - (neg_f[i] / sum_Q);
  }
  free(pos_f);
  free(neg_f);
  delete tree;
}


// Evaluate t-SNE cost function (approximately)
double TSNE::evaluateError(size_t* row_P, size_t* col_P,
  double* val_P, double* Y, size_t N, size_t D, double theta) {
  // Get estimate of normalization term
  SPTree* tree = new SPTree(D, Y, N);
  static std::unique_ptr<double[]> buff (new double[D]);
  double sum_Q = .0;
  for (size_t n = 0; n < N; n++) {
    tree->computeNonEdgeForces(n, theta, buff.get(), &sum_Q);
  }

  // Loop over all edges to compute t-SNE error
  size_t ind1, ind2;
  double C = .0, Q;
  #pragma omp parallel for schedule(dynamic)
  for (size_t n = 0; n < N; n++) {
    ind1 = n * D;
    for (size_t i = row_P[n]; i < row_P[n + 1]; i++) {
      Q = .0;
      ind2 = col_P[i] * D;
      for (size_t d = 0; d < D; d++) buff[d]  = Y[ind1 + d];
      for (size_t d = 0; d < D; d++) buff[d] -= Y[ind2 + d];
      for (size_t d = 0; d < D; d++) Q += buff[d] * buff[d];
      Q = (1.0 / (1.0 + Q)) / sum_Q;
      C += val_P[i] * log((val_P[i] + FLT_MIN) / (Q + FLT_MIN));
    }
  }
  // Clean up memory
  delete tree;
  return C;
}


// Compute input similarities with a fixed perplexity using ball trees
// (this function allocates memory another function should free)
void TSNE::computeGaussianPerplexity(double* X, size_t N, size_t D, size_t* row_P,
  size_t* col_P, double* val_P, double perplexity, size_t K) {
    if (perplexity > static_cast<double>(K)) {
      printf("Perplexity should be lower than K!\n");
    }
    std::unique_ptr<double[]>cur_P (new double[N-1]);
    if (cur_P == NULL) {
      printf("Memory allocation failed!\n");
      exit(1);
    }
    row_P[0] = 0;
    for (size_t n = 0; n < N; n++) {
      row_P[n + 1] = row_P[n] + static_cast<size_t> (K);
    }

    // Build ball tree on data set
    VpTree<DataPoint, DataPoint::euclidean_distance>* tree =
      new VpTree<DataPoint, DataPoint::euclidean_distance>();
    std::vector<DataPoint> obj_X(N, DataPoint(D, -1, X));

    #pragma omp parallel for schedule(dynamic)
    for (size_t n = 0; n < N; n++) {
      obj_X[n] = DataPoint(D, static_cast<int>(n), X + n * D);
    }
    tree->create(obj_X);

    // Loop over all points to find nearest neighbors
    printf("Building tree...\n");
    std::vector<DataPoint> indices;
    std::vector<double> distances;

    #pragma omp parallel for schedule(dynamic) private(indices, distances)
    for (size_t n = 0; n < N; n++) {
      if (n % 10000 == 0) {
        printf(" - point %zu of %zu\n", n, N);
      }
      // Find nearest neighbors
      indices.clear();
      distances.clear();
      tree->search(obj_X[n], K + 1, &indices, &distances);

      // Initialize some variables for binary search
      bool found = false;
      double beta = 1.0;
      double min_beta = -DBL_MAX;
      double max_beta =  DBL_MAX;
      double tol = 1e-5;

      // Iterate until we found a good perplexity
      int iter = 0; double sum_P;
      while (!found && iter < 200) {
        // Compute Gaussian kernel row
        for (size_t m = 0; m < K; m++) {
          cur_P[m] = exp(-beta * distances[m + 1] * distances[m + 1]);
        }

        // Compute entropy of current row
        sum_P = DBL_MIN;
        for (size_t m = 0; m < K; m++) {
          sum_P += cur_P[m];
        }
        double H = .0;
        for (size_t m = 0; m < K; m++) {
          H += beta * (distances[m + 1] * distances[m + 1] * cur_P[m]);
        }
        H = (H / sum_P) + log(sum_P);

        // Evaluate whether the entropy is within the tolerance level
        double Hdiff = H - log(perplexity);
        if (Hdiff < tol && -Hdiff < tol) {
          found = true;
        } else {
          if (Hdiff > 0) {
            min_beta = beta;
            if (max_beta == DBL_MAX || max_beta == -DBL_MAX) {
                beta *= 2.0;
            } else {
              beta = (beta + max_beta) / 2.0;
            }
          } else {
            max_beta = beta;
            if (min_beta == -DBL_MAX || min_beta == DBL_MAX) {
                beta /= 2.0;
            } else {
              beta = (beta + min_beta) / 2.0;
            }
          }
        }
        // Update iteration counter
        iter++;
      }

      // Row-normalize current row of P and store in matrix
      #pragma omp parallel for schedule(dynamic)
      for (size_t m = 0; m < K; m++) {
        cur_P[m] /= sum_P;
        col_P[row_P[n] + m] = (size_t) indices[m + 1].index();
        val_P[row_P[n] + m] = cur_P[m];
      }
    }
    // Clean up memory
    obj_X.clear();
    delete tree;
}


// Symmetrizes a sparse matrix
void TSNE::symmetrizeMatrix(std::unique_ptr<size_t[]> &_row_P,
  std::unique_ptr<size_t[]> &_col_P,
  std::unique_ptr<double[]> &_val_P, size_t N) {
  // Get sparse matrix
  size_t* row_P = _row_P.get();
  size_t* col_P = _col_P.get();
  double* val_P = _val_P.get();

  // Count number of elements and row counts of symmetric matrix
  int* row_counts = reinterpret_cast<int*> (calloc(N, sizeof(int)));
  if (row_counts == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }

  #pragma omp parallel for schedule(dynamic)
  for (size_t n = 0; n < N; n++) {
    for (size_t i = row_P[n]; i < row_P[n + 1]; i++) {
      // Check whether element (col_P[i], n) is present
      bool present = false;
      for (size_t m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
        if (col_P[m] == n) present = true;
      }
      if (present) {
        row_counts[n]++;
      } else {
        row_counts[n]++;
        row_counts[col_P[i]]++;
      }
    }
  }
  int no_elem = 0;
  for (size_t n = 0; n < N; n++) {
    no_elem += row_counts[n];
  }

  // Allocate memory for symmetrized matrix
  size_t* sym_row_P = (size_t*) malloc((N + 1) * sizeof(size_t));
  size_t* sym_col_P = (size_t*) malloc(no_elem * sizeof(size_t));
  double* sym_val_P = reinterpret_cast<double*> (malloc(no_elem * sizeof(double)));

  if (sym_row_P == NULL || sym_col_P == NULL || sym_val_P == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }

  // Construct new row indices for symmetric matrix
  sym_row_P[0] = 0;
  for (size_t n = 0; n < N; n++) {
    sym_row_P[n + 1] = sym_row_P[n] + (size_t) row_counts[n];
  }

  // Fill the result matrix
  int* offset = reinterpret_cast<int*> (calloc(N, sizeof(int)));

  for (size_t n = 0; n < N; n++) {
    for (size_t i = row_P[n]; i < row_P[n + 1]; i++) {  // considering element(n, col_P[i])
      // Check whether element (col_P[i], n) is present
      bool present = false;
      for (size_t m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
        if (col_P[m] == n) {
          present = true;
          if (n <= col_P[i]) {  // make sure we do not add elements twice
            sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
            sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
            sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i] + val_P[m];
            sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
          }
        }
      }

      // If (col_P[i], n) is not present, there is no addition involved
      if (!present) {
        sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
        sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
        sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i];
        sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
      }

      // Update offsets
      if (!present || (present && n <= col_P[i])) {
        offset[n]++;
        if (col_P[i] != n) offset[col_P[i]]++;
      }
    }
  }

  // Divide the result by two
  for (int i = 0; i < no_elem; i++) {
    sym_val_P[i] /= 2.0;
  }

  // Return symmetrized matrices
  _row_P.reset(sym_row_P);
  _col_P.reset(sym_col_P);
  _val_P.reset(sym_val_P);

  // Free up some memory
  free(offset); offset = NULL;
  free(row_counts); row_counts  = NULL;
}

// Compute squared Euclidean distance matrix
void TSNE::computeSquaredEuclideanDistance(double* X, size_t N, size_t D, double* DD) {
  const double* XnD = X;
  for (size_t n = 0; n < N; ++n, XnD += D) {
    const double* XmD = XnD + D;
    double* curr_elem = &DD[n*N + n];
    *curr_elem = 0.0;
    double* curr_elem_sym = curr_elem + N;
    for (size_t m = n + 1; m < N; ++m, XmD+=D, curr_elem_sym+=N) {
      *(++curr_elem) = 0.0;
      for (size_t d = 0; d < D; ++d) {
          *curr_elem += (XnD[d] - XmD[d]) * (XnD[d] - XmD[d]);
      }
      *curr_elem_sym = *curr_elem;
    }
  }
}


// Makes data zero-mean
void TSNE::zeroMean(double* X, size_t N, size_t D) {
  // Compute data mean
 static std::unique_ptr<double[]> mean (new double[D]);
  size_t nD = 0;
  #pragma omp parallel for schedule(dynamic)
  for (size_t n = 0; n < N; n++) {
    for (size_t d = 0; d < D; d++) {
      mean[d] += X[nD + d];
    }
    nD += D;
  }

  #pragma omp parallel for schedule(dynamic)
  for (size_t d = 0; d < D; d++) {
    mean[d] /= static_cast<double> (N);
  }

  // Subtract data mean
  nD = 0;
  for (size_t n = 0; n < N; n++) {
    for (size_t d = 0; d < D; d++) {
      X[nD + d] -= mean[d];
    }
    nD += D;
  }
}


// Generates a Gaussian random number
double TSNE::randn() {
  double x, y, radius;
  do {
    x = 2 * (rand() / (static_cast<double> (RAND_MAX) + 1)) - 1;
    y = 2 * (rand() / (static_cast<double> (RAND_MAX) + 1)) - 1;
    radius = (x * x) + (y * y);
  } while ((radius >= 1.0) || (radius == 0.0));
  radius = sqrt(-2 * log(radius) / radius);
  x *= radius;
  y *= radius;
  return x;
}

}  // namespace datadriven
}  // namespace sgpp
