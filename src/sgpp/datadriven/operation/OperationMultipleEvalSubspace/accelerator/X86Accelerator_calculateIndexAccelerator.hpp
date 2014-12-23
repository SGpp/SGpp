static inline void calculateIndexAccelerator(size_t dim, size_t nextIterationToRecalc,
					  const double * const (&dataTuplePtr)[8], 
					  uint32_t *hInversePtr,
					  uint32_t *(&intermediates)[8],
					  double *(&evalIndexValues)[8], 
					  uint32_t (&indexFlat)[8],
					  double (&phiEval)[8]
					  ) {

  indexFlat[0] = intermediates[0][nextIterationToRecalc];
  indexFlat[1] = intermediates[1][nextIterationToRecalc];
  indexFlat[2] = intermediates[2][nextIterationToRecalc];
  indexFlat[3] = intermediates[3][nextIterationToRecalc];

  indexFlat[4] = intermediates[4][nextIterationToRecalc];
  indexFlat[5] = intermediates[5][nextIterationToRecalc];
  indexFlat[6] = intermediates[6][nextIterationToRecalc];
  indexFlat[7] = intermediates[7][nextIterationToRecalc];

  phiEval[0] = evalIndexValues[0][nextIterationToRecalc];
  phiEval[1] = evalIndexValues[1][nextIterationToRecalc];
  phiEval[2] = evalIndexValues[2][nextIterationToRecalc];
  phiEval[3] = evalIndexValues[3][nextIterationToRecalc];

  phiEval[4] = evalIndexValues[4][nextIterationToRecalc];
  phiEval[5] = evalIndexValues[5][nextIterationToRecalc];
  phiEval[6] = evalIndexValues[6][nextIterationToRecalc];
  phiEval[7] = evalIndexValues[7][nextIterationToRecalc];

  uint32_t mask = 0x1;

  for (size_t i = nextIterationToRecalc; i < dim; i++) {
    double unadjusted[8];
    unadjusted[0] = dataTuplePtr[0][i] * hInversePtr[i];
    unadjusted[1] = dataTuplePtr[1][i] * hInversePtr[i];
    unadjusted[2] = dataTuplePtr[2][i] * hInversePtr[i];
    unadjusted[3] = dataTuplePtr[3][i] * hInversePtr[i];

    unadjusted[4] = dataTuplePtr[4][i] * hInversePtr[i];
    unadjusted[5] = dataTuplePtr[5][i] * hInversePtr[i];
    unadjusted[6] = dataTuplePtr[6][i] * hInversePtr[i];
    unadjusted[7] = dataTuplePtr[7][i] * hInversePtr[i];
	  
    //implies flooring
    uint32_t rounded[8];
    rounded[0] = static_cast<uint32_t>(unadjusted[0]);
    rounded[1] = static_cast<uint32_t>(unadjusted[1]);
    rounded[2] = static_cast<uint32_t>(unadjusted[2]);
    rounded[3] = static_cast<uint32_t>(unadjusted[3]);

    rounded[4] = static_cast<uint32_t>(unadjusted[4]);
    rounded[5] = static_cast<uint32_t>(unadjusted[5]);
    rounded[6] = static_cast<uint32_t>(unadjusted[6]);
    rounded[7] = static_cast<uint32_t>(unadjusted[7]);
	  
    uint32_t sign[8];
    sign[0] = mask ^ (mask & rounded[0]);
    sign[1] = mask ^ (mask & rounded[1]);
    sign[2] = mask ^ (mask & rounded[2]);
    sign[3] = mask ^ (mask & rounded[3]);

    sign[4] = mask ^ (mask & rounded[4]);
    sign[5] = mask ^ (mask & rounded[5]);
    sign[6] = mask ^ (mask & rounded[6]);
    sign[7] = mask ^ (mask & rounded[7]);

    uint32_t index[8];
    index[0] = rounded[0] + sign[0];
    index[1] = rounded[1] + sign[1];
    index[2] = rounded[2] + sign[2];
    index[3] = rounded[3] + sign[3];

    index[4] = rounded[4] + sign[4];
    index[5] = rounded[5] + sign[5];
    index[6] = rounded[6] + sign[6];
    index[7] = rounded[7] + sign[7];

    // indexPtr[0][i] = rounded[0] + sign[0];
    // indexPtr[1][i] = rounded[1] + sign[1];
    // indexPtr[2][i] = rounded[2] + sign[2];
    // indexPtr[3][i] = rounded[3] + sign[3];
    
    //flatten index
    int actualDirectionGridPoints = hInversePtr[i];
    actualDirectionGridPoints >>= 1;

    indexFlat[0] *= actualDirectionGridPoints;
    indexFlat[1] *= actualDirectionGridPoints;
    indexFlat[2] *= actualDirectionGridPoints;
    indexFlat[3] *= actualDirectionGridPoints;

    indexFlat[4] *= actualDirectionGridPoints;
    indexFlat[5] *= actualDirectionGridPoints;
    indexFlat[6] *= actualDirectionGridPoints;
    indexFlat[7] *= actualDirectionGridPoints;

    size_t actualIndex[8];
    // actualIndex[0] = indexPtr[0][i];
    // actualIndex[1] = indexPtr[1][i];
    // actualIndex[2] = indexPtr[2][i];
    // actualIndex[3] = indexPtr[3][i];

    actualIndex[0] = index[0];
    actualIndex[1] = index[1];
    actualIndex[2] = index[2];
    actualIndex[3] = index[3];

    actualIndex[4] = index[4];
    actualIndex[5] = index[5];
    actualIndex[6] = index[6];
    actualIndex[7] = index[7];

    actualIndex[0] >>= 1; //divide index by 2, skip even indices
    actualIndex[1] >>= 1;
    actualIndex[2] >>= 1;
    actualIndex[3] >>= 1;

    actualIndex[4] >>= 1;
    actualIndex[5] >>= 1;
    actualIndex[6] >>= 1;
    actualIndex[7] >>= 1;

    indexFlat[0] += actualIndex[0];
    indexFlat[1] += actualIndex[1];
    indexFlat[2] += actualIndex[2];
    indexFlat[3] += actualIndex[3];

    indexFlat[4] += actualIndex[4];
    indexFlat[5] += actualIndex[5];
    indexFlat[6] += actualIndex[6];
    indexFlat[7] += actualIndex[7];

    intermediates[0][i + 1] = indexFlat[0];
    intermediates[1][i + 1] = indexFlat[1];
    intermediates[2][i + 1] = indexFlat[2];
    intermediates[3][i + 1] = indexFlat[3];

    intermediates[4][i + 1] = indexFlat[4];
    intermediates[5][i + 1] = indexFlat[5];
    intermediates[6][i + 1] = indexFlat[6];
    intermediates[7][i + 1] = indexFlat[7];
    
    //eval basis function
    double phi1DEval[8];
    phi1DEval[0] = hInversePtr[i] * dataTuplePtr[0][i] - index[0];
    phi1DEval[1] = hInversePtr[i] * dataTuplePtr[1][i] - index[1];
    phi1DEval[2] = hInversePtr[i] * dataTuplePtr[2][i] - index[2];
    phi1DEval[3] = hInversePtr[i] * dataTuplePtr[3][i] - index[3];

    phi1DEval[4] = hInversePtr[i] * dataTuplePtr[4][i] - index[4];
    phi1DEval[5] = hInversePtr[i] * dataTuplePtr[5][i] - index[5];
    phi1DEval[6] = hInversePtr[i] * dataTuplePtr[6][i] - index[6];
    phi1DEval[7] = hInversePtr[i] * dataTuplePtr[7][i] - index[7];

    phi1DEval[0] = max(0.0, 1.0 - abs(phi1DEval[0]));
    phi1DEval[1] = max(0.0, 1.0 - abs(phi1DEval[1]));
    phi1DEval[2] = max(0.0, 1.0 - abs(phi1DEval[2]));
    phi1DEval[3] = max(0.0, 1.0 - abs(phi1DEval[3]));

    phi1DEval[4] = max(0.0, 1.0 - abs(phi1DEval[4]));
    phi1DEval[5] = max(0.0, 1.0 - abs(phi1DEval[5]));
    phi1DEval[6] = max(0.0, 1.0 - abs(phi1DEval[6]));
    phi1DEval[7] = max(0.0, 1.0 - abs(phi1DEval[7]));

    phiEval[0] *= phi1DEval[0];
    phiEval[1] *= phi1DEval[1];
    phiEval[2] *= phi1DEval[2];
    phiEval[3] *= phi1DEval[3];

    phiEval[4] *= phi1DEval[4];
    phiEval[5] *= phi1DEval[5];
    phiEval[6] *= phi1DEval[6];
    phiEval[7] *= phi1DEval[7];

    evalIndexValues[0][i + 1] = phiEval[0];
    evalIndexValues[1][i + 1] = phiEval[1];
    evalIndexValues[2][i + 1] = phiEval[2];
    evalIndexValues[3][i + 1] = phiEval[3];

    evalIndexValues[4][i + 1] = phiEval[4];
    evalIndexValues[5][i + 1] = phiEval[5];
    evalIndexValues[6][i + 1] = phiEval[6];
    evalIndexValues[7][i + 1] = phiEval[7];
  }
}

//TODO: required until transpose is also ported to accelerator data structure
static inline void calculateIndexAccelerator(size_t dim, size_t nextIterationToRecalc,
                      const double * const (&dataTuplePtr)[8],
                      vector<uint32_t> &hInversePtr,
                      uint32_t *(&intermediates)[8],
                      double *(&evalIndexValues)[8],
                      uint32_t (&indexFlat)[8],
                      double (&phiEval)[8]
                      ) {

  indexFlat[0] = intermediates[0][nextIterationToRecalc];
  indexFlat[1] = intermediates[1][nextIterationToRecalc];
  indexFlat[2] = intermediates[2][nextIterationToRecalc];
  indexFlat[3] = intermediates[3][nextIterationToRecalc];

  indexFlat[4] = intermediates[4][nextIterationToRecalc];
  indexFlat[5] = intermediates[5][nextIterationToRecalc];
  indexFlat[6] = intermediates[6][nextIterationToRecalc];
  indexFlat[7] = intermediates[7][nextIterationToRecalc];

  phiEval[0] = evalIndexValues[0][nextIterationToRecalc];
  phiEval[1] = evalIndexValues[1][nextIterationToRecalc];
  phiEval[2] = evalIndexValues[2][nextIterationToRecalc];
  phiEval[3] = evalIndexValues[3][nextIterationToRecalc];

  phiEval[4] = evalIndexValues[4][nextIterationToRecalc];
  phiEval[5] = evalIndexValues[5][nextIterationToRecalc];
  phiEval[6] = evalIndexValues[6][nextIterationToRecalc];
  phiEval[7] = evalIndexValues[7][nextIterationToRecalc];

  uint32_t mask = 0x1;

  for (size_t i = nextIterationToRecalc; i < dim; i++) {
    double unadjusted[8];
    unadjusted[0] = dataTuplePtr[0][i] * hInversePtr[i];
    unadjusted[1] = dataTuplePtr[1][i] * hInversePtr[i];
    unadjusted[2] = dataTuplePtr[2][i] * hInversePtr[i];
    unadjusted[3] = dataTuplePtr[3][i] * hInversePtr[i];

    unadjusted[4] = dataTuplePtr[4][i] * hInversePtr[i];
    unadjusted[5] = dataTuplePtr[5][i] * hInversePtr[i];
    unadjusted[6] = dataTuplePtr[6][i] * hInversePtr[i];
    unadjusted[7] = dataTuplePtr[7][i] * hInversePtr[i];

    //implies flooring
    uint32_t rounded[8];
    rounded[0] = static_cast<uint32_t>(unadjusted[0]);
    rounded[1] = static_cast<uint32_t>(unadjusted[1]);
    rounded[2] = static_cast<uint32_t>(unadjusted[2]);
    rounded[3] = static_cast<uint32_t>(unadjusted[3]);

    rounded[4] = static_cast<uint32_t>(unadjusted[4]);
    rounded[5] = static_cast<uint32_t>(unadjusted[5]);
    rounded[6] = static_cast<uint32_t>(unadjusted[6]);
    rounded[7] = static_cast<uint32_t>(unadjusted[7]);

    uint32_t sign[8];
    sign[0] = mask ^ (mask & rounded[0]);
    sign[1] = mask ^ (mask & rounded[1]);
    sign[2] = mask ^ (mask & rounded[2]);
    sign[3] = mask ^ (mask & rounded[3]);

    sign[4] = mask ^ (mask & rounded[4]);
    sign[5] = mask ^ (mask & rounded[5]);
    sign[6] = mask ^ (mask & rounded[6]);
    sign[7] = mask ^ (mask & rounded[7]);

    uint32_t index[8];
    index[0] = rounded[0] + sign[0];
    index[1] = rounded[1] + sign[1];
    index[2] = rounded[2] + sign[2];
    index[3] = rounded[3] + sign[3];

    index[4] = rounded[4] + sign[4];
    index[5] = rounded[5] + sign[5];
    index[6] = rounded[6] + sign[6];
    index[7] = rounded[7] + sign[7];

    // indexPtr[0][i] = rounded[0] + sign[0];
    // indexPtr[1][i] = rounded[1] + sign[1];
    // indexPtr[2][i] = rounded[2] + sign[2];
    // indexPtr[3][i] = rounded[3] + sign[3];

    //flatten index
    int actualDirectionGridPoints = hInversePtr[i];
    actualDirectionGridPoints >>= 1;

    indexFlat[0] *= actualDirectionGridPoints;
    indexFlat[1] *= actualDirectionGridPoints;
    indexFlat[2] *= actualDirectionGridPoints;
    indexFlat[3] *= actualDirectionGridPoints;

    indexFlat[4] *= actualDirectionGridPoints;
    indexFlat[5] *= actualDirectionGridPoints;
    indexFlat[6] *= actualDirectionGridPoints;
    indexFlat[7] *= actualDirectionGridPoints;

    size_t actualIndex[8];
    // actualIndex[0] = indexPtr[0][i];
    // actualIndex[1] = indexPtr[1][i];
    // actualIndex[2] = indexPtr[2][i];
    // actualIndex[3] = indexPtr[3][i];

    actualIndex[0] = index[0];
    actualIndex[1] = index[1];
    actualIndex[2] = index[2];
    actualIndex[3] = index[3];

    actualIndex[4] = index[4];
    actualIndex[5] = index[5];
    actualIndex[6] = index[6];
    actualIndex[7] = index[7];

    actualIndex[0] >>= 1; //divide index by 2, skip even indices
    actualIndex[1] >>= 1;
    actualIndex[2] >>= 1;
    actualIndex[3] >>= 1;

    actualIndex[4] >>= 1;
    actualIndex[5] >>= 1;
    actualIndex[6] >>= 1;
    actualIndex[7] >>= 1;

    indexFlat[0] += actualIndex[0];
    indexFlat[1] += actualIndex[1];
    indexFlat[2] += actualIndex[2];
    indexFlat[3] += actualIndex[3];

    indexFlat[4] += actualIndex[4];
    indexFlat[5] += actualIndex[5];
    indexFlat[6] += actualIndex[6];
    indexFlat[7] += actualIndex[7];

    intermediates[0][i + 1] = indexFlat[0];
    intermediates[1][i + 1] = indexFlat[1];
    intermediates[2][i + 1] = indexFlat[2];
    intermediates[3][i + 1] = indexFlat[3];

    intermediates[4][i + 1] = indexFlat[4];
    intermediates[5][i + 1] = indexFlat[5];
    intermediates[6][i + 1] = indexFlat[6];
    intermediates[7][i + 1] = indexFlat[7];

    //eval basis function
    double phi1DEval[8];
    phi1DEval[0] = hInversePtr[i] * dataTuplePtr[0][i] - index[0];
    phi1DEval[1] = hInversePtr[i] * dataTuplePtr[1][i] - index[1];
    phi1DEval[2] = hInversePtr[i] * dataTuplePtr[2][i] - index[2];
    phi1DEval[3] = hInversePtr[i] * dataTuplePtr[3][i] - index[3];

    phi1DEval[4] = hInversePtr[i] * dataTuplePtr[4][i] - index[4];
    phi1DEval[5] = hInversePtr[i] * dataTuplePtr[5][i] - index[5];
    phi1DEval[6] = hInversePtr[i] * dataTuplePtr[6][i] - index[6];
    phi1DEval[7] = hInversePtr[i] * dataTuplePtr[7][i] - index[7];

    phi1DEval[0] = max(0.0, 1.0 - abs(phi1DEval[0]));
    phi1DEval[1] = max(0.0, 1.0 - abs(phi1DEval[1]));
    phi1DEval[2] = max(0.0, 1.0 - abs(phi1DEval[2]));
    phi1DEval[3] = max(0.0, 1.0 - abs(phi1DEval[3]));

    phi1DEval[4] = max(0.0, 1.0 - abs(phi1DEval[4]));
    phi1DEval[5] = max(0.0, 1.0 - abs(phi1DEval[5]));
    phi1DEval[6] = max(0.0, 1.0 - abs(phi1DEval[6]));
    phi1DEval[7] = max(0.0, 1.0 - abs(phi1DEval[7]));

    phiEval[0] *= phi1DEval[0];
    phiEval[1] *= phi1DEval[1];
    phiEval[2] *= phi1DEval[2];
    phiEval[3] *= phi1DEval[3];

    phiEval[4] *= phi1DEval[4];
    phiEval[5] *= phi1DEval[5];
    phiEval[6] *= phi1DEval[6];
    phiEval[7] *= phi1DEval[7];

    evalIndexValues[0][i + 1] = phiEval[0];
    evalIndexValues[1][i + 1] = phiEval[1];
    evalIndexValues[2][i + 1] = phiEval[2];
    evalIndexValues[3][i + 1] = phiEval[3];

    evalIndexValues[4][i + 1] = phiEval[4];
    evalIndexValues[5][i + 1] = phiEval[5];
    evalIndexValues[6][i + 1] = phiEval[6];
    evalIndexValues[7][i + 1] = phiEval[7];
  }
}

// static inline void calculateIndexAccelerator(size_t dim, size_t nextIterationToRecalc,
// 					  const double * const (&dataTuplePtr)[4], 
// 					  vector<uint32_t> &hInversePtr, 
// 					  uint32_t *(&intermediates)[4],
// 					  double *(&evalIndexValues)[4], 
// 					  //uint32_t *(&indexPtr)[4],
// 					  uint32_t (&indexFlat)[4],
// 					  double (&phiEval)[4]
// 					  ) {

//   indexFlat[0] = intermediates[0][nextIterationToRecalc];
//   indexFlat[1] = intermediates[1][nextIterationToRecalc];
//   indexFlat[2] = intermediates[2][nextIterationToRecalc];
//   indexFlat[3] = intermediates[3][nextIterationToRecalc];

//   phiEval[0] = evalIndexValues[0][nextIterationToRecalc];
//   phiEval[1] = evalIndexValues[1][nextIterationToRecalc];
//   phiEval[2] = evalIndexValues[2][nextIterationToRecalc];
//   phiEval[3] = evalIndexValues[3][nextIterationToRecalc];

//   uint32_t mask = 0x1;

//   for (size_t i = nextIterationToRecalc; i < dim; i++) {
//     double unadjusted[4];
//     unadjusted[0] = dataTuplePtr[0][i] * hInversePtr[i];
//     unadjusted[1] = dataTuplePtr[1][i] * hInversePtr[i];
//     unadjusted[2] = dataTuplePtr[2][i] * hInversePtr[i];
//     unadjusted[3] = dataTuplePtr[3][i] * hInversePtr[i];
	  
//     //implies flooring
//     uint32_t rounded[4];
//     rounded[0] = static_cast<uint32_t>(unadjusted[0]);
//     rounded[1] = static_cast<uint32_t>(unadjusted[1]);
//     rounded[2] = static_cast<uint32_t>(unadjusted[2]);
//     rounded[3] = static_cast<uint32_t>(unadjusted[3]);	 
	  
//     uint32_t sign[4];
//     sign[0] = mask ^ (mask & rounded[0]);
//     sign[1] = mask ^ (mask & rounded[1]);
//     sign[2] = mask ^ (mask & rounded[2]);
//     sign[3] = mask ^ (mask & rounded[3]);

//     uint32_t index[4];
//     index[0] = rounded[0] + sign[0];
//     index[1] = rounded[1] + sign[1];
//     index[2] = rounded[2] + sign[2];
//     index[3] = rounded[3] + sign[3];      

//     // indexPtr[0][i] = rounded[0] + sign[0];
//     // indexPtr[1][i] = rounded[1] + sign[1];
//     // indexPtr[2][i] = rounded[2] + sign[2];
//     // indexPtr[3][i] = rounded[3] + sign[3];
    
//     //flatten index
//     int actualDirectionGridPoints = hInversePtr[i];
//     actualDirectionGridPoints >>= 1;

//     indexFlat[0] *= actualDirectionGridPoints;
//     indexFlat[1] *= actualDirectionGridPoints;
//     indexFlat[2] *= actualDirectionGridPoints;
//     indexFlat[3] *= actualDirectionGridPoints;

//     size_t actualIndex[4];
//     // actualIndex[0] = indexPtr[0][i];
//     // actualIndex[1] = indexPtr[1][i];
//     // actualIndex[2] = indexPtr[2][i];
//     // actualIndex[3] = indexPtr[3][i];

//     actualIndex[0] = index[0];
//     actualIndex[1] = index[1];
//     actualIndex[2] = index[2];
//     actualIndex[3] = index[3];


//     actualIndex[0] >>= 1; //divide index by 2, skip even indices
//     actualIndex[1] >>= 1;
//     actualIndex[2] >>= 1;
//     actualIndex[3] >>= 1;

//     indexFlat[0] += actualIndex[0];
//     indexFlat[1] += actualIndex[1];
//     indexFlat[2] += actualIndex[2];
//     indexFlat[3] += actualIndex[3];

//     intermediates[0][i + 1] = indexFlat[0];
//     intermediates[1][i + 1] = indexFlat[1];
//     intermediates[2][i + 1] = indexFlat[2];
//     intermediates[3][i + 1] = indexFlat[3];
    
//     //eval basis function
//     double phi1DEval[4];
//     phi1DEval[0] = hInversePtr[i] * dataTuplePtr[0][i] - index[0];
//     phi1DEval[1] = hInversePtr[i] * dataTuplePtr[1][i] - index[1];
//     phi1DEval[2] = hInversePtr[i] * dataTuplePtr[2][i] - index[2];
//     phi1DEval[3] = hInversePtr[i] * dataTuplePtr[3][i] - index[3];

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

// static inline void calculateIndexAccelerator(size_t dim, size_t nextIterationToRecalc,
// 					  const double * const (&dataTuplePtr)[4], 
// 					  vector<uint32_t> &hInversePtr, 
// 					  uint32_t *(&intermediates)[4],
// 					  double *(&evalIndexValues)[4], 
// 					  //uint32_t *(&indexPtr)[4],
// 					  uint32_t (&indexFlat)[4],
// 					  double (&phiEval)[4]
// 					  ) {

//   __m128i oneIntegerReg = _mm_set1_epi32((uint32_t) 1);

//   union {
//     __m128d doubleRegister;
//     __m128i integerRegister;
//     uint32_t uint32Value[4];
//   } sseUnion;

//   // flatten only
//   __m128i indexFlatReg = _mm_set_epi32(intermediates[3][nextIterationToRecalc], 
// 				       intermediates[2][nextIterationToRecalc], 
// 				       intermediates[1][nextIterationToRecalc], 
// 				       intermediates[0][nextIterationToRecalc]);


//   // evaluate only
//   union {
//     __m256d doubleRegister;
//     double doubleValue[4];
//   } avxUnion;

//   int64_t absIMask = 0x7FFFFFFFFFFFFFFF;
//   double* fabsMask = (double *) &absIMask;
//   __m256d absMask = _mm256_broadcast_sd(fabsMask);
//   __m256d one = _mm256_set1_pd(1.0);
//   __m256d zero = _mm256_set1_pd(0.0);

//   __m256d phiEvalReg = _mm256_set_pd(evalIndexValues[3][nextIterationToRecalc], 
//   				     evalIndexValues[2][nextIterationToRecalc], 
//   				     evalIndexValues[1][nextIterationToRecalc], 
//   				     evalIndexValues[0][nextIterationToRecalc]);

//   for (size_t i = nextIterationToRecalc; i < dim; i += 1) {


//   //for (size_t outer = nextIterationToRecalc; outer < dim; outer += 4) {
    
//     // __m256d dataReg0 = _mm256_loadu_pd(dataTuplePtr[0] + outer);
//     // __m256d dataReg1 = _mm256_loadu_pd(dataTuplePtr[1] + outer);
//     // __m256d dataReg2 = _mm256_loadu_pd(dataTuplePtr[2] + outer);
//     // __m256d dataReg3 = _mm256_loadu_pd(dataTuplePtr[3] + outer);
//     // __m256d t0low = _mm256_unpacklo_pd(dataReg0, dataReg1);
//     // __m256d t0hi = _mm256_unpackhi_pd(dataReg0, dataReg1);
//     // __m256d t1low = _mm256_unpacklo_pd(dataReg2, dataReg3);
//     // __m256d t1hi = _mm256_unpackhi_pd(dataReg2, dataReg3);

//     // __m256d allData[4];
//     // allData[0] = _mm256_permute2f128_pd(t0low, t1low, 0 + 32);
//     // allData[1] = _mm256_permute2f128_pd(t0hi, t1hi, 0 + 32);
//     // allData[2] = _mm256_permute2f128_pd(t0low, t1low, 1 + 48);
//     // allData[3] = _mm256_permute2f128_pd(t0hi, t1hi, 1 + 48);

//     // size_t dataVectorIndex = 0;

//     //for (size_t i = outer; i < min(outer + 4, dim); i++) {
//       //TODO: replace by avx2 gather
//       //TODO (better): replace by array access
//       __m256d dataTupleReg = _mm256_set_pd(dataTuplePtr[3][i], dataTuplePtr[2][i], 
//        					   dataTuplePtr[1][i], dataTuplePtr[0][i]);      

//       // __m256d dataTupleReg = allData[dataVectorIndex];
//       // dataVectorIndex += 1;
           
//       //TODO: replace by array access
//       __m256d hInverseReg = _mm256_set1_pd((double) hInversePtr[i]);
//       __m256d unadjustedReg = _mm256_mul_pd(dataTupleReg, hInverseReg);

//       //implies flooring
//       __m128i roundedReg = _mm256_cvttpd_epi32(unadjustedReg);	 
//       __m128i andedReg = _mm_and_si128(oneIntegerReg, roundedReg);
//       __m128i signReg = _mm_xor_si128(oneIntegerReg, andedReg);
//       __m128i indexReg = _mm_add_epi32(roundedReg, signReg);	  	  

//       //flatten index
//       //TODO: can be precomputed as array
//       uint32_t actualDirectionGridPoints = hInversePtr[i];
//       actualDirectionGridPoints >>= 1;
//       __m128i actualDirectionGridPointsReg = _mm_set1_epi32(actualDirectionGridPoints);

//       indexFlatReg = _mm_mullo_epi32(indexFlatReg, actualDirectionGridPointsReg);

//       __m128i indexShiftedReg = _mm_srli_epi32(indexReg, 1);
	  
//       indexFlatReg = _mm_add_epi32(indexFlatReg, indexShiftedReg);
    
//       sseUnion.integerRegister = indexFlatReg;
//       intermediates[0][i + 1] = sseUnion.uint32Value[0];
//       intermediates[1][i + 1] = sseUnion.uint32Value[1];
//       intermediates[2][i + 1] = sseUnion.uint32Value[2];
//       intermediates[3][i + 1] = sseUnion.uint32Value[3];

//       // evaluate
//       __m256d indexDoubleReg = _mm256_cvtepi32_pd(indexReg);

//       __m256d phi1DEvalReg = _mm256_mul_pd(hInverseReg, dataTupleReg);
//       phi1DEvalReg = _mm256_sub_pd(phi1DEvalReg, indexDoubleReg);

//       phi1DEvalReg = _mm256_and_pd(phi1DEvalReg, absMask);
//       phi1DEvalReg = _mm256_sub_pd(one, phi1DEvalReg);
//       phi1DEvalReg = _mm256_max_pd(zero, phi1DEvalReg);

//       phiEvalReg = _mm256_mul_pd(phiEvalReg, phi1DEvalReg);
    
//       avxUnion.doubleRegister = phiEvalReg;
//       evalIndexValues[0][i + 1] = avxUnion.doubleValue[0];
//       evalIndexValues[1][i + 1] = avxUnion.doubleValue[1];
//       evalIndexValues[2][i + 1] = avxUnion.doubleValue[2];
//       evalIndexValues[3][i + 1] = avxUnion.doubleValue[3];
//   }
//   //may a structure ind[0] im[0] eval[0] null ind[1] ... might help
//   //}
//   _mm_storeu_si128((__m128i *) indexFlat, indexFlatReg);
//   _mm256_storeu_pd(phiEval, phiEvalReg);
// }
