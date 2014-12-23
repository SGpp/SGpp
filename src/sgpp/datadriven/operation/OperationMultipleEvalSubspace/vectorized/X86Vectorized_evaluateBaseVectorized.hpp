// static inline void evaluateBaseVectorized(size_t dim, size_t nextIterationToRecalc,
// 					  double *(&dataTuplePtr)[4], 
// 					  double *(&evalIndexValues)[4], 
// 					  uint32_t *hInversePtr, uint32_t *(&indexPtr)[4],
// 					  double (&phiEval)[4]) {
//   phiEval[0] = evalIndexValues[0][nextIterationToRecalc];
//   phiEval[1] = evalIndexValues[1][nextIterationToRecalc];
//   phiEval[2] = evalIndexValues[2][nextIterationToRecalc];
//   phiEval[3] = evalIndexValues[3][nextIterationToRecalc];

//   //prepare the values for the individual components
//   for (size_t i = nextIterationToRecalc; i < dim; i++) {
//     double phi1DEval[4];
//     phi1DEval[0] = hInversePtr[i] * dataTuplePtr[0][i] - indexPtr[0][i];
//     phi1DEval[1] = hInversePtr[i] * dataTuplePtr[1][i] - indexPtr[1][i];
//     phi1DEval[2] = hInversePtr[i] * dataTuplePtr[2][i] - indexPtr[2][i];
//     phi1DEval[3] = hInversePtr[i] * dataTuplePtr[3][i] - indexPtr[3][i];

//     phi1DEval[0] = max(0.0, 1.0 - abs(phi1DEval[0]));
//     phi1DEval[1] = max(0.0, 1.0 - abs(phi1DEval[1]));
//     phi1DEval[2] = max(0.0, 1.0 - abs(phi1DEval[2]));
//     phi1DEval[3] = max(0.0, 1.0 - abs(phi1DEval[3]));

//     phiEval[0] *= phi1DEval[0];
//     phiEval[1] *= phi1DEval[1];
//     phiEval[2] *= phi1DEval[2];
//     phiEval[3] *= phi1DEval[3];

//     evalIndexValues[0][i + 1] = phiEval[0];
//     evalIndexValues[1][i + 1] = phiEval[1];
//     evalIndexValues[2][i + 1] = phiEval[2];
//     evalIndexValues[3][i + 1] = phiEval[3];
//   }
// }

static inline void evaluateBaseVectorized(size_t dim, size_t nextIterationToRecalc,
					  double *(&dataTuplePtr)[4], 
					  double *(&evalIndexValues)[4], 
					  uint32_t *hInversePtr, uint32_t *(&indexPtr)[4],
					  double (&phiEval)[4]) {
  union {
    __m256d doubleRegister;
    double doubleValue[4];
  } avxUnion;

  int64_t absIMask = 0x7FFFFFFFFFFFFFFF;
  double* fabsMask = (double *) &absIMask;
  __m256d absMask = _mm256_broadcast_sd(fabsMask);
  __m256d one = _mm256_set1_pd(1.0);
  __m256d zero = _mm256_set1_pd(0.0);

  __m256d phiEvalReg = _mm256_set_pd(evalIndexValues[3][nextIterationToRecalc], 
				     evalIndexValues[2][nextIterationToRecalc], 
				     evalIndexValues[1][nextIterationToRecalc], 
				     evalIndexValues[0][nextIterationToRecalc]);

  //prepare the values for the individual components
  for (size_t i = nextIterationToRecalc; i < dim; i++) {
    __m256d hInverseReg = _mm256_set1_pd((double) hInversePtr[i]);
    
    __m256d dataTupleReg = _mm256_set_pd(dataTuplePtr[3][i], dataTuplePtr[2][i], 
					 dataTuplePtr[1][i], dataTuplePtr[0][i]);

    __m128i indexIntegerReg = _mm_set_epi32(indexPtr[3][i], indexPtr[2][i],
				     indexPtr[1][i], indexPtr[0][i]);

    __m256d indexReg = _mm256_cvtepi32_pd(indexIntegerReg);

    __m256d phi1DEvalReg = _mm256_mul_pd(hInverseReg, dataTupleReg);
    phi1DEvalReg = _mm256_sub_pd(phi1DEvalReg, indexReg);

    phi1DEvalReg = _mm256_and_pd(phi1DEvalReg, absMask);
    phi1DEvalReg = _mm256_sub_pd(one, phi1DEvalReg);
    phi1DEvalReg = _mm256_max_pd(zero, phi1DEvalReg);

    phiEvalReg = _mm256_mul_pd(phiEvalReg, phi1DEvalReg);
    
    avxUnion.doubleRegister = phiEvalReg;
    evalIndexValues[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues[3][i + 1] = avxUnion.doubleValue[3];
  }

  _mm256_storeu_pd(phiEval, phiEvalReg);
}
