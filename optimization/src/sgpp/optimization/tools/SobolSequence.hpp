// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// This Code was copied from https://people.sc.fsu.edu/~jburkardt/cpp_src/sobol/sobol.html
// (visited on 28.11.2018). It is distributed under the GNU LGPL License

#include <string>

namespace sgpp {
namespace optimization {

int i4_bit_hi1(int n);
int i4_bit_lo0(int n);
int i4_max(int i1, int i2);
int i4_min(int i1, int i2);
void i4_sobol(int dim_num, int *seed, float quasi[]);
float *i4_sobol_generate(int m, int n, int skip);
int i4_uniform(int b, int c, int *seed);

int i8_bit_hi1(long long int n);
int i8_bit_lo0(long long int n);
long long int i8_max(long long int i1, long long int i2);
long long int i8_min(long long int i1, long long int i2);
/**
 * return the point of the sobol sequence according to the seed and number of dimensiosn
 */
void i8_sobol(int dim_num, long long int *seed, double quasi[]);
double *i8_sobol_generate(int m, int n, int skip);
long long int i8_uniform(long long int b, long long int c, int *seed);

float r4_abs(float x);
int r4_nint(float x);
float r4_uniform_01(int *seed);

double r8_abs(double x);
int r8_nint(double x);
double r8_uniform_01(int *seed);

void r8mat_write(std::string output_filename, int m, int n, double table[]);

/**
 * returns a seed which apparently is a good fit for the given number of dimensions
 */
int tau_sobol(int dim_num);
void sobol_timestamp();
}  // namespace optimization
}  // namespace sgpp
