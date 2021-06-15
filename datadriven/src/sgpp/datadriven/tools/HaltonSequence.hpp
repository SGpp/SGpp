// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// This Code was copied from https://people.sc.fsu.edu/~jburkardt/cpp_src/halton/halton.html
// (visited on 28.11.2018). It is distributed under the GNU LGPL License

#include <string>

namespace sgpp {
namespace datadriven {

double *halton(int i, int m);
/**
 * gives the i'th point in the m-dimensional halton squence according to the array of primes b
 */
double *halton_base(int i, int m, int b[]);
int halton_inverse(double r[], int m);
double *halton_sequence(int i1, int i2, int m);
int i4vec_sum(int n, int a[]);
/**
 * return the n'th prime number. used to create the array of primes for halton_base
 */
int prime(int n);
double r8_mod(double x, double y);
void r8mat_print(int m, int n, double a[], std::string title);
void r8mat_print_some(int m, int n, double a[], int ilo, int jlo, int ihi, int jhi,
                      std::string title);
void halton_timestamp();
}  // namespace datadriven
}  // namespace sgpp
