// static inline void calculateIndexVectorized(size_t dim, size_t nextIterationToRecalc, 
// 					    double *(&dataTuplePtr)[4], 
// 					    uint32_t *hInversePtr, uint32_t *(&indexPtr)[4]) {
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
	  
//     uint32_t mask = 0x1;
	  
//     uint32_t sign[4];
//     sign[0] = mask ^ (mask & rounded[0]);
//     sign[1] = mask ^ (mask & rounded[1]);
//     sign[2] = mask ^ (mask & rounded[2]);
//     sign[3] = mask ^ (mask & rounded[3]);
	  
//     indexPtr[0][i] = rounded[0] + sign[0];
//     indexPtr[1][i] = rounded[1] + sign[1];
//     indexPtr[2][i] = rounded[2] + sign[2];
//     indexPtr[3][i] = rounded[3] + sign[3];
//   }
// }

static inline void calculateIndexVectorized(size_t dim, size_t nextIterationToRecalc, 
					    const double * const (&dataTuplePtr)[4], 
					    uint32_t *hInversePtr, 
					    uint32_t *(&intermediates)[4],
					    double *(&evalIndexValues)[4], 
					    //uint32_t *(&indexPtr)[4],
					    uint32_t (&indexFlat)[4],
					    double (&phiEval)[4]
					    ) {

  __m128i oneIntegerReg = _mm_set1_epi32((uint32_t) 1);

  union {
    __m128d doubleRegister;
    __m128i integerRegister;
    uint32_t uint32Value[4];
  } sseUnion;

  // flatten only
  __m128i indexFlatReg = _mm_set_epi32(intermediates[3][nextIterationToRecalc], 
				       intermediates[2][nextIterationToRecalc], 
				       intermediates[1][nextIterationToRecalc], 
				       intermediates[0][nextIterationToRecalc]);

  // evaluate only
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

  for (size_t i = nextIterationToRecalc; i < dim; i += 1) {
  //for (size_t outer = nextIterationToRecalc; outer < dim; outer += 4) {
    
    // __m256d dataReg0 = _mm256_loadu_pd(dataTuplePtr[0] + outer);
    // __m256d dataReg1 = _mm256_loadu_pd(dataTuplePtr[1] + outer);
    // __m256d dataReg2 = _mm256_loadu_pd(dataTuplePtr[2] + outer);
    // __m256d dataReg3 = _mm256_loadu_pd(dataTuplePtr[3] + outer);
    // __m256d t0low = _mm256_unpacklo_pd(dataReg0, dataReg1);
    // __m256d t0hi = _mm256_unpackhi_pd(dataReg0, dataReg1);
    // __m256d t1low = _mm256_unpacklo_pd(dataReg2, dataReg3);
    // __m256d t1hi = _mm256_unpackhi_pd(dataReg2, dataReg3);

    // __m256d allData[4];
    // allData[0] = _mm256_permute2f128_pd(t0low, t1low, 0 + 32);
    // allData[1] = _mm256_permute2f128_pd(t0hi, t1hi, 0 + 32);
    // allData[2] = _mm256_permute2f128_pd(t0low, t1low, 1 + 48);
    // allData[3] = _mm256_permute2f128_pd(t0hi, t1hi, 1 + 48);

    // size_t dataVectorIndex = 0;

    //for (size_t i = outer; i < min(outer + 4, dim); i++) {
      //TODO: replace by avx2 gather
      //TODO (better): replace by array access
      __m256d dataTupleReg = _mm256_set_pd(dataTuplePtr[3][i], dataTuplePtr[2][i], 
       					   dataTuplePtr[1][i], dataTuplePtr[0][i]);      

      // __m256d dataTupleReg = allData[dataVectorIndex];
      // dataVectorIndex += 1;
           
      //TODO: replace by array access
      __m256d hInverseReg = _mm256_set1_pd((double) hInversePtr[i]);
      __m256d unadjustedReg = _mm256_mul_pd(dataTupleReg, hInverseReg);

      //implies flooring
      __m128i roundedReg = _mm256_cvttpd_epi32(unadjustedReg);	 
      __m128i andedReg = _mm_and_si128(oneIntegerReg, roundedReg);
      __m128i signReg = _mm_xor_si128(oneIntegerReg, andedReg);
      __m128i indexReg = _mm_add_epi32(roundedReg, signReg);	  	  

      //flatten index
      //TODO: can be precomputed as array
      uint32_t actualDirectionGridPoints = hInversePtr[i];
      actualDirectionGridPoints >>= 1;
      __m128i actualDirectionGridPointsReg = _mm_set1_epi32(actualDirectionGridPoints);

      indexFlatReg = _mm_mullo_epi32(indexFlatReg, actualDirectionGridPointsReg);

      __m128i indexShiftedReg = _mm_srli_epi32(indexReg, 1);
	  
      indexFlatReg = _mm_add_epi32(indexFlatReg, indexShiftedReg);
    
      sseUnion.integerRegister = indexFlatReg;
      intermediates[0][i + 1] = sseUnion.uint32Value[0];
      intermediates[1][i + 1] = sseUnion.uint32Value[1];
      intermediates[2][i + 1] = sseUnion.uint32Value[2];
      intermediates[3][i + 1] = sseUnion.uint32Value[3];

      // evaluate
      __m256d indexDoubleReg = _mm256_cvtepi32_pd(indexReg);

      __m256d phi1DEvalReg = _mm256_mul_pd(hInverseReg, dataTupleReg);
      phi1DEvalReg = _mm256_sub_pd(phi1DEvalReg, indexDoubleReg);

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
  //may a structure ind[0] im[0] eval[0] null ind[1] ... might help
  //}
  _mm_storeu_si128((__m128i *) indexFlat, indexFlatReg);
  _mm256_storeu_pd(phiEval, phiEvalReg);
}

static inline void calculateIndexVectorized2(size_t dim, size_t nextIterationToRecalc, 
					     //rep
					     const double * const (&dataTuplePtr)[4], 
					     const double * const (&dataTuplePtr2)[4], 
					     uint32_t *hInversePtr, 
					     //rep
					     uint32_t *(&intermediates)[4],
					     uint32_t *(&intermediates2)[4],
					     //rep
					     double *(&evalIndexValues)[4],
					     double *(&evalIndexValues2)[4],
					     //rep
					     uint32_t (&indexFlat)[4],
					     uint32_t (&indexFlat2)[4],
					     //rep
					     double (&phiEval)[4],
					     double (&phiEval2)[4]) {

  __m128i oneIntegerReg = _mm_set1_epi32((uint32_t) 1);

  union {
    __m128d doubleRegister;
    __m128i integerRegister;
    uint32_t uint32Value[4];
  } sseUnion;

  // flatten only
  __m128i indexFlatReg = _mm_set_epi32(intermediates[3][nextIterationToRecalc], 
				       intermediates[2][nextIterationToRecalc], 
				       intermediates[1][nextIterationToRecalc], 
				       intermediates[0][nextIterationToRecalc]);
  __m128i indexFlatReg2 = _mm_set_epi32(intermediates2[3][nextIterationToRecalc], 
					intermediates2[2][nextIterationToRecalc], 
					intermediates2[1][nextIterationToRecalc], 
					intermediates2[0][nextIterationToRecalc]);

  // evaluate only
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
  __m256d phiEvalReg2 = _mm256_set_pd(evalIndexValues2[3][nextIterationToRecalc], 
				      evalIndexValues2[2][nextIterationToRecalc], 
				      evalIndexValues2[1][nextIterationToRecalc], 
				      evalIndexValues2[0][nextIterationToRecalc]);

  for (size_t i = nextIterationToRecalc; i < dim; i += 1) {
    //TODO: replace by avx2 gather
    //TODO (better): replace by array access
    __m256d dataTupleReg = _mm256_set_pd(dataTuplePtr[3][i], dataTuplePtr[2][i], 
					 dataTuplePtr[1][i], dataTuplePtr[0][i]);      
    __m256d dataTupleReg2 = _mm256_set_pd(dataTuplePtr2[3][i], dataTuplePtr2[2][i], 
					  dataTuplePtr2[1][i], dataTuplePtr2[0][i]);      

    // __m256d dataTupleReg = allData[dataVectorIndex];
    // dataVectorIndex += 1;
           
    //TODO: replace by array access
    __m256d hInverseReg = _mm256_set1_pd((double) hInversePtr[i]);

    __m256d unadjustedReg = _mm256_mul_pd(dataTupleReg, hInverseReg);
    __m256d unadjustedReg2 = _mm256_mul_pd(dataTupleReg2, hInverseReg);

    //implies flooring
    __m128i roundedReg = _mm256_cvttpd_epi32(unadjustedReg);	 
    __m128i roundedReg2 = _mm256_cvttpd_epi32(unadjustedReg2);
	 
    __m128i andedReg = _mm_and_si128(oneIntegerReg, roundedReg);
    __m128i andedReg2 = _mm_and_si128(oneIntegerReg, roundedReg2);

    __m128i signReg = _mm_xor_si128(oneIntegerReg, andedReg);
    __m128i signReg2 = _mm_xor_si128(oneIntegerReg, andedReg2);

    __m128i indexReg = _mm_add_epi32(roundedReg, signReg);	  	  
    __m128i indexReg2 = _mm_add_epi32(roundedReg2, signReg2);

    //flatten index
    //TODO: can be precomputed as array
    uint32_t actualDirectionGridPoints = hInversePtr[i];
    actualDirectionGridPoints >>= 1;
    __m128i actualDirectionGridPointsReg = _mm_set1_epi32(actualDirectionGridPoints);

    indexFlatReg = _mm_mullo_epi32(indexFlatReg, actualDirectionGridPointsReg);
    indexFlatReg2 = _mm_mullo_epi32(indexFlatReg2, actualDirectionGridPointsReg);

    __m128i indexShiftedReg = _mm_srli_epi32(indexReg, 1);
    __m128i indexShiftedReg2 = _mm_srli_epi32(indexReg2, 1);
	  
    indexFlatReg = _mm_add_epi32(indexFlatReg, indexShiftedReg);
    indexFlatReg2 = _mm_add_epi32(indexFlatReg2, indexShiftedReg2);
    
    sseUnion.integerRegister = indexFlatReg;
    intermediates[0][i + 1] = sseUnion.uint32Value[0];
    intermediates[1][i + 1] = sseUnion.uint32Value[1];
    intermediates[2][i + 1] = sseUnion.uint32Value[2];
    intermediates[3][i + 1] = sseUnion.uint32Value[3];

    sseUnion.integerRegister = indexFlatReg2;
    intermediates2[0][i + 1] = sseUnion.uint32Value[0];
    intermediates2[1][i + 1] = sseUnion.uint32Value[1];
    intermediates2[2][i + 1] = sseUnion.uint32Value[2];
    intermediates2[3][i + 1] = sseUnion.uint32Value[3];

    // evaluate
    __m256d indexDoubleReg = _mm256_cvtepi32_pd(indexReg);
    __m256d indexDoubleReg2 = _mm256_cvtepi32_pd(indexReg2);

    __m256d phi1DEvalReg = _mm256_mul_pd(hInverseReg, dataTupleReg);
    __m256d phi1DEvalReg2 = _mm256_mul_pd(hInverseReg, dataTupleReg2);

    phi1DEvalReg = _mm256_sub_pd(phi1DEvalReg, indexDoubleReg);
    phi1DEvalReg2 = _mm256_sub_pd(phi1DEvalReg2, indexDoubleReg2);

    phi1DEvalReg = _mm256_and_pd(phi1DEvalReg, absMask);
    phi1DEvalReg2 = _mm256_and_pd(phi1DEvalReg2, absMask);

    phi1DEvalReg = _mm256_sub_pd(one, phi1DEvalReg);
    phi1DEvalReg2 = _mm256_sub_pd(one, phi1DEvalReg2);

    phi1DEvalReg = _mm256_max_pd(zero, phi1DEvalReg);
    phi1DEvalReg2 = _mm256_max_pd(zero, phi1DEvalReg2);

    phiEvalReg = _mm256_mul_pd(phiEvalReg, phi1DEvalReg);
    phiEvalReg2 = _mm256_mul_pd(phiEvalReg2, phi1DEvalReg2);
    
    avxUnion.doubleRegister = phiEvalReg;
    evalIndexValues[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues[3][i + 1] = avxUnion.doubleValue[3];

    avxUnion.doubleRegister = phiEvalReg2;
    evalIndexValues2[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues2[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues2[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues2[3][i + 1] = avxUnion.doubleValue[3];
  }
  //may a structure ind[0] im[0] eval[0] null ind[1] ... might help
  _mm_storeu_si128((__m128i *) indexFlat, indexFlatReg);
  _mm_storeu_si128((__m128i *) indexFlat2, indexFlatReg2);
    
  _mm256_storeu_pd(phiEval, phiEvalReg);
  _mm256_storeu_pd(phiEval2, phiEvalReg2);
}

static inline void calculateIndexVectorized6(size_t dim, size_t nextIterationToRecalc, 
					     //rep
					     const double * const (&dataTuplePtr)[4], 
					     const double * const (&dataTuplePtr2)[4], 
					     const double * const (&dataTuplePtr3)[4], 
					     const double * const (&dataTuplePtr4)[4], 
					     const double * const (&dataTuplePtr5)[4], 
					     const double * const (&dataTuplePtr6)[4], 
					     uint32_t *hInversePtr, 
					     //rep
					     uint32_t *(&intermediates)[4],
					     uint32_t *(&intermediates2)[4],
					     uint32_t *(&intermediates3)[4],
					     uint32_t *(&intermediates4)[4],
					     uint32_t *(&intermediates5)[4],
					     uint32_t *(&intermediates6)[4],
					     //rep
					     double *(&evalIndexValues)[4],
					     double *(&evalIndexValues2)[4],
					     double *(&evalIndexValues3)[4],
					     double *(&evalIndexValues4)[4],
					     double *(&evalIndexValues5)[4],
					     double *(&evalIndexValues6)[4],
					     //rep
					     uint32_t (&indexFlat)[4],
					     uint32_t (&indexFlat2)[4],
					     uint32_t (&indexFlat3)[4],
					     uint32_t (&indexFlat4)[4],
					     uint32_t (&indexFlat5)[4],
					     uint32_t (&indexFlat6)[4],
					     //rep
					     double (&phiEval)[4],
					     double (&phiEval2)[4],
					     double (&phiEval3)[4],
					     double (&phiEval4)[4],
					     double (&phiEval5)[4],
					     double (&phiEval6)[4]
					     ) {

  __m128i oneIntegerReg = _mm_set1_epi32((uint32_t) 1);

  union {
    __m128d doubleRegister;
    __m128i integerRegister;
    uint32_t uint32Value[4];
  } sseUnion;

  // flatten only
  __m128i indexFlatReg = _mm_set_epi32(intermediates[3][nextIterationToRecalc], 
				       intermediates[2][nextIterationToRecalc], 
				       intermediates[1][nextIterationToRecalc], 
				       intermediates[0][nextIterationToRecalc]);
  __m128i indexFlatReg2 = _mm_set_epi32(intermediates2[3][nextIterationToRecalc], 
					intermediates2[2][nextIterationToRecalc], 
					intermediates2[1][nextIterationToRecalc], 
					intermediates2[0][nextIterationToRecalc]);
  __m128i indexFlatReg3 = _mm_set_epi32(intermediates3[3][nextIterationToRecalc], 
					intermediates3[2][nextIterationToRecalc], 
					intermediates3[1][nextIterationToRecalc], 
					intermediates3[0][nextIterationToRecalc]);
  __m128i indexFlatReg4 = _mm_set_epi32(intermediates4[3][nextIterationToRecalc], 
					intermediates4[2][nextIterationToRecalc], 
					intermediates4[1][nextIterationToRecalc], 
					intermediates4[0][nextIterationToRecalc]);
  __m128i indexFlatReg5 = _mm_set_epi32(intermediates5[3][nextIterationToRecalc], 
					intermediates5[2][nextIterationToRecalc], 
					intermediates5[1][nextIterationToRecalc], 
					intermediates5[0][nextIterationToRecalc]);
  __m128i indexFlatReg6 = _mm_set_epi32(intermediates6[3][nextIterationToRecalc], 
					intermediates6[2][nextIterationToRecalc], 
					intermediates6[1][nextIterationToRecalc], 
					intermediates6[0][nextIterationToRecalc]);

  // evaluate only
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
  __m256d phiEvalReg2 = _mm256_set_pd(evalIndexValues2[3][nextIterationToRecalc], 
				      evalIndexValues2[2][nextIterationToRecalc], 
				      evalIndexValues2[1][nextIterationToRecalc], 
				      evalIndexValues2[0][nextIterationToRecalc]);
  __m256d phiEvalReg3 = _mm256_set_pd(evalIndexValues3[3][nextIterationToRecalc], 
				      evalIndexValues3[2][nextIterationToRecalc], 
				      evalIndexValues3[1][nextIterationToRecalc], 
				      evalIndexValues3[0][nextIterationToRecalc]);
  __m256d phiEvalReg4 = _mm256_set_pd(evalIndexValues4[3][nextIterationToRecalc], 
				      evalIndexValues4[2][nextIterationToRecalc], 
				      evalIndexValues4[1][nextIterationToRecalc], 
				      evalIndexValues4[0][nextIterationToRecalc]);
  __m256d phiEvalReg5 = _mm256_set_pd(evalIndexValues5[3][nextIterationToRecalc], 
				      evalIndexValues5[2][nextIterationToRecalc], 
				      evalIndexValues5[1][nextIterationToRecalc], 
				      evalIndexValues5[0][nextIterationToRecalc]);
  __m256d phiEvalReg6 = _mm256_set_pd(evalIndexValues6[3][nextIterationToRecalc], 
				      evalIndexValues6[2][nextIterationToRecalc], 
				      evalIndexValues6[1][nextIterationToRecalc], 
				      evalIndexValues6[0][nextIterationToRecalc]);

  for (size_t i = nextIterationToRecalc; i < dim; i += 1) {
    //TODO: replace by avx2 gather
    //TODO (better): replace by array access
    __m256d dataTupleReg = _mm256_set_pd(dataTuplePtr[3][i], dataTuplePtr[2][i], 
					 dataTuplePtr[1][i], dataTuplePtr[0][i]);      
    __m256d dataTupleReg2 = _mm256_set_pd(dataTuplePtr2[3][i], dataTuplePtr2[2][i], 
					  dataTuplePtr2[1][i], dataTuplePtr2[0][i]);
    __m256d dataTupleReg3 = _mm256_set_pd(dataTuplePtr3[3][i], dataTuplePtr3[2][i], 
					  dataTuplePtr3[1][i], dataTuplePtr3[0][i]);
    __m256d dataTupleReg4 = _mm256_set_pd(dataTuplePtr4[3][i], dataTuplePtr4[2][i], 
					  dataTuplePtr4[1][i], dataTuplePtr4[0][i]);
    __m256d dataTupleReg5 = _mm256_set_pd(dataTuplePtr5[3][i], dataTuplePtr5[2][i], 
					  dataTuplePtr5[1][i], dataTuplePtr5[0][i]);
    __m256d dataTupleReg6 = _mm256_set_pd(dataTuplePtr6[3][i], dataTuplePtr6[2][i], 
					  dataTuplePtr6[1][i], dataTuplePtr6[0][i]);

    // __m256d dataTupleReg = allData[dataVectorIndex];
    // dataVectorIndex += 1;
           
    //TODO: replace by array access
    __m256d hInverseReg = _mm256_set1_pd((double) hInversePtr[i]);

    __m256d unadjustedReg = _mm256_mul_pd(dataTupleReg, hInverseReg);
    __m256d unadjustedReg2 = _mm256_mul_pd(dataTupleReg2, hInverseReg);
    __m256d unadjustedReg3 = _mm256_mul_pd(dataTupleReg3, hInverseReg);
    __m256d unadjustedReg4 = _mm256_mul_pd(dataTupleReg4, hInverseReg);
    __m256d unadjustedReg5 = _mm256_mul_pd(dataTupleReg5, hInverseReg);
    __m256d unadjustedReg6 = _mm256_mul_pd(dataTupleReg6, hInverseReg);

    //implies flooring
    __m128i roundedReg = _mm256_cvttpd_epi32(unadjustedReg);	 
    __m128i roundedReg2 = _mm256_cvttpd_epi32(unadjustedReg2);
    __m128i roundedReg3 = _mm256_cvttpd_epi32(unadjustedReg3);
    __m128i roundedReg4 = _mm256_cvttpd_epi32(unadjustedReg4);
    __m128i roundedReg5 = _mm256_cvttpd_epi32(unadjustedReg5);
    __m128i roundedReg6 = _mm256_cvttpd_epi32(unadjustedReg6);
	 
    __m128i andedReg = _mm_and_si128(oneIntegerReg, roundedReg);
    __m128i andedReg2 = _mm_and_si128(oneIntegerReg, roundedReg2);
    __m128i andedReg3 = _mm_and_si128(oneIntegerReg, roundedReg3);
    __m128i andedReg4 = _mm_and_si128(oneIntegerReg, roundedReg4);
    __m128i andedReg5 = _mm_and_si128(oneIntegerReg, roundedReg5);
    __m128i andedReg6 = _mm_and_si128(oneIntegerReg, roundedReg6);

    __m128i signReg = _mm_xor_si128(oneIntegerReg, andedReg);
    __m128i signReg2 = _mm_xor_si128(oneIntegerReg, andedReg2);
    __m128i signReg3 = _mm_xor_si128(oneIntegerReg, andedReg3);
    __m128i signReg4 = _mm_xor_si128(oneIntegerReg, andedReg4);
    __m128i signReg5 = _mm_xor_si128(oneIntegerReg, andedReg5);
    __m128i signReg6 = _mm_xor_si128(oneIntegerReg, andedReg6);

    __m128i indexReg = _mm_add_epi32(roundedReg, signReg);	  	  
    __m128i indexReg2 = _mm_add_epi32(roundedReg2, signReg2);
    __m128i indexReg3 = _mm_add_epi32(roundedReg3, signReg3);
    __m128i indexReg4 = _mm_add_epi32(roundedReg4, signReg4);
    __m128i indexReg5 = _mm_add_epi32(roundedReg5, signReg5);
    __m128i indexReg6 = _mm_add_epi32(roundedReg6, signReg6);

    //flatten index
    //TODO: can be precomputed as array
    uint32_t actualDirectionGridPoints = hInversePtr[i];
    actualDirectionGridPoints >>= 1;
    __m128i actualDirectionGridPointsReg = _mm_set1_epi32(actualDirectionGridPoints);

    indexFlatReg = _mm_mullo_epi32(indexFlatReg, actualDirectionGridPointsReg);
    indexFlatReg2 = _mm_mullo_epi32(indexFlatReg2, actualDirectionGridPointsReg);
    indexFlatReg3 = _mm_mullo_epi32(indexFlatReg3, actualDirectionGridPointsReg);
    indexFlatReg4 = _mm_mullo_epi32(indexFlatReg4, actualDirectionGridPointsReg);
    indexFlatReg5 = _mm_mullo_epi32(indexFlatReg5, actualDirectionGridPointsReg);
    indexFlatReg6 = _mm_mullo_epi32(indexFlatReg6, actualDirectionGridPointsReg);

    __m128i indexShiftedReg = _mm_srli_epi32(indexReg, 1);
    __m128i indexShiftedReg2 = _mm_srli_epi32(indexReg2, 1);
    __m128i indexShiftedReg3 = _mm_srli_epi32(indexReg3, 1);
    __m128i indexShiftedReg4 = _mm_srli_epi32(indexReg4, 1);
    __m128i indexShiftedReg5 = _mm_srli_epi32(indexReg5, 1);
    __m128i indexShiftedReg6 = _mm_srli_epi32(indexReg6, 1);
	  
    indexFlatReg = _mm_add_epi32(indexFlatReg, indexShiftedReg);
    indexFlatReg2 = _mm_add_epi32(indexFlatReg2, indexShiftedReg2);
    indexFlatReg3 = _mm_add_epi32(indexFlatReg3, indexShiftedReg3);
    indexFlatReg4 = _mm_add_epi32(indexFlatReg4, indexShiftedReg4);
    indexFlatReg5 = _mm_add_epi32(indexFlatReg5, indexShiftedReg5);
    indexFlatReg6 = _mm_add_epi32(indexFlatReg6, indexShiftedReg6);
    
    sseUnion.integerRegister = indexFlatReg;
    intermediates[0][i + 1] = sseUnion.uint32Value[0];
    intermediates[1][i + 1] = sseUnion.uint32Value[1];
    intermediates[2][i + 1] = sseUnion.uint32Value[2];
    intermediates[3][i + 1] = sseUnion.uint32Value[3];

    sseUnion.integerRegister = indexFlatReg2;
    intermediates2[0][i + 1] = sseUnion.uint32Value[0];
    intermediates2[1][i + 1] = sseUnion.uint32Value[1];
    intermediates2[2][i + 1] = sseUnion.uint32Value[2];
    intermediates2[3][i + 1] = sseUnion.uint32Value[3];

    sseUnion.integerRegister = indexFlatReg3;
    intermediates3[0][i + 1] = sseUnion.uint32Value[0];
    intermediates3[1][i + 1] = sseUnion.uint32Value[1];
    intermediates3[2][i + 1] = sseUnion.uint32Value[2];
    intermediates3[3][i + 1] = sseUnion.uint32Value[3];

    sseUnion.integerRegister = indexFlatReg4;
    intermediates4[0][i + 1] = sseUnion.uint32Value[0];
    intermediates4[1][i + 1] = sseUnion.uint32Value[1];
    intermediates4[2][i + 1] = sseUnion.uint32Value[2];
    intermediates4[3][i + 1] = sseUnion.uint32Value[3];

    sseUnion.integerRegister = indexFlatReg5;
    intermediates5[0][i + 1] = sseUnion.uint32Value[0];
    intermediates5[1][i + 1] = sseUnion.uint32Value[1];
    intermediates5[2][i + 1] = sseUnion.uint32Value[2];
    intermediates5[3][i + 1] = sseUnion.uint32Value[3];

    sseUnion.integerRegister = indexFlatReg6;
    intermediates6[0][i + 1] = sseUnion.uint32Value[0];
    intermediates6[1][i + 1] = sseUnion.uint32Value[1];
    intermediates6[2][i + 1] = sseUnion.uint32Value[2];
    intermediates6[3][i + 1] = sseUnion.uint32Value[3];

    // evaluate
    __m256d indexDoubleReg = _mm256_cvtepi32_pd(indexReg);
    __m256d indexDoubleReg2 = _mm256_cvtepi32_pd(indexReg2);
    __m256d indexDoubleReg3 = _mm256_cvtepi32_pd(indexReg3);
    __m256d indexDoubleReg4 = _mm256_cvtepi32_pd(indexReg4);
    __m256d indexDoubleReg5 = _mm256_cvtepi32_pd(indexReg5);
    __m256d indexDoubleReg6 = _mm256_cvtepi32_pd(indexReg6);

    __m256d phi1DEvalReg = _mm256_mul_pd(hInverseReg, dataTupleReg);
    __m256d phi1DEvalReg2 = _mm256_mul_pd(hInverseReg, dataTupleReg2);
    __m256d phi1DEvalReg3 = _mm256_mul_pd(hInverseReg, dataTupleReg3);
    __m256d phi1DEvalReg4 = _mm256_mul_pd(hInverseReg, dataTupleReg4);
    __m256d phi1DEvalReg5 = _mm256_mul_pd(hInverseReg, dataTupleReg5);
    __m256d phi1DEvalReg6 = _mm256_mul_pd(hInverseReg, dataTupleReg6);

    phi1DEvalReg = _mm256_sub_pd(phi1DEvalReg, indexDoubleReg);
    phi1DEvalReg2 = _mm256_sub_pd(phi1DEvalReg2, indexDoubleReg2);
    phi1DEvalReg3 = _mm256_sub_pd(phi1DEvalReg3, indexDoubleReg3);
    phi1DEvalReg4 = _mm256_sub_pd(phi1DEvalReg4, indexDoubleReg4);
    phi1DEvalReg5 = _mm256_sub_pd(phi1DEvalReg5, indexDoubleReg5);
    phi1DEvalReg6 = _mm256_sub_pd(phi1DEvalReg6, indexDoubleReg6);

    phi1DEvalReg = _mm256_and_pd(phi1DEvalReg, absMask);
    phi1DEvalReg2 = _mm256_and_pd(phi1DEvalReg2, absMask);
    phi1DEvalReg3 = _mm256_and_pd(phi1DEvalReg3, absMask);
    phi1DEvalReg4 = _mm256_and_pd(phi1DEvalReg4, absMask);
    phi1DEvalReg5 = _mm256_and_pd(phi1DEvalReg5, absMask);
    phi1DEvalReg6 = _mm256_and_pd(phi1DEvalReg6, absMask);

    phi1DEvalReg = _mm256_sub_pd(one, phi1DEvalReg);
    phi1DEvalReg2 = _mm256_sub_pd(one, phi1DEvalReg2);
    phi1DEvalReg3 = _mm256_sub_pd(one, phi1DEvalReg3);
    phi1DEvalReg4 = _mm256_sub_pd(one, phi1DEvalReg4);
    phi1DEvalReg5 = _mm256_sub_pd(one, phi1DEvalReg5);
    phi1DEvalReg6 = _mm256_sub_pd(one, phi1DEvalReg6);

    phi1DEvalReg = _mm256_max_pd(zero, phi1DEvalReg);
    phi1DEvalReg2 = _mm256_max_pd(zero, phi1DEvalReg2);
    phi1DEvalReg3 = _mm256_max_pd(zero, phi1DEvalReg3);
    phi1DEvalReg4 = _mm256_max_pd(zero, phi1DEvalReg4);
    phi1DEvalReg5 = _mm256_max_pd(zero, phi1DEvalReg5);
    phi1DEvalReg6 = _mm256_max_pd(zero, phi1DEvalReg6);

    phiEvalReg = _mm256_mul_pd(phiEvalReg, phi1DEvalReg);
    phiEvalReg2 = _mm256_mul_pd(phiEvalReg2, phi1DEvalReg2);
    phiEvalReg3 = _mm256_mul_pd(phiEvalReg3, phi1DEvalReg3);
    phiEvalReg4 = _mm256_mul_pd(phiEvalReg4, phi1DEvalReg4);
    phiEvalReg5 = _mm256_mul_pd(phiEvalReg5, phi1DEvalReg5);
    phiEvalReg6 = _mm256_mul_pd(phiEvalReg6, phi1DEvalReg6);
    
    avxUnion.doubleRegister = phiEvalReg;
    evalIndexValues[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues[3][i + 1] = avxUnion.doubleValue[3];

    avxUnion.doubleRegister = phiEvalReg2;
    evalIndexValues2[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues2[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues2[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues2[3][i + 1] = avxUnion.doubleValue[3];

    avxUnion.doubleRegister = phiEvalReg3;
    evalIndexValues3[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues3[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues3[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues3[3][i + 1] = avxUnion.doubleValue[3];

    avxUnion.doubleRegister = phiEvalReg4;
    evalIndexValues4[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues4[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues4[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues4[3][i + 1] = avxUnion.doubleValue[3];

    avxUnion.doubleRegister = phiEvalReg5;
    evalIndexValues5[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues5[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues5[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues5[3][i + 1] = avxUnion.doubleValue[3];

    avxUnion.doubleRegister = phiEvalReg6;
    evalIndexValues6[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues6[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues6[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues6[3][i + 1] = avxUnion.doubleValue[3];
  }
  //may a structure ind[0] im[0] eval[0] null ind[1] ... might help
  _mm_storeu_si128((__m128i *) indexFlat, indexFlatReg);
  _mm_storeu_si128((__m128i *) indexFlat2, indexFlatReg2);
  _mm_storeu_si128((__m128i *) indexFlat3, indexFlatReg3);
  _mm_storeu_si128((__m128i *) indexFlat4, indexFlatReg4);
  _mm_storeu_si128((__m128i *) indexFlat5, indexFlatReg5);
  _mm_storeu_si128((__m128i *) indexFlat6, indexFlatReg6);
    
  _mm256_storeu_pd(phiEval, phiEvalReg);
  _mm256_storeu_pd(phiEval2, phiEvalReg2);
  _mm256_storeu_pd(phiEval3, phiEvalReg3);
  _mm256_storeu_pd(phiEval4, phiEvalReg4);
  _mm256_storeu_pd(phiEval5, phiEvalReg5);
  _mm256_storeu_pd(phiEval6, phiEvalReg6);
}

/*
static inline void calculateIndexVectorized(size_t dim, size_t nextIterationToRecalc, 
					    double *(&dataTuplePtr)[4], uint32_t *hInversePtr, 
					    uint32_t *(&intermediates)[4],
					    double *(&evalIndexValues)[4], 
					    //uint32_t *(&indexPtr)[4],
					    uint32_t (&indexFlat)[4],
					    double (&phiEval)[4]
					    ) {

  __m128i oneIntegerReg = _mm_set1_epi32((uint32_t) 1);

  union {
    __m128d doubleRegister;
    __m128i integerRegister;
    uint32_t uint32Value[4];
  } sseUnion;

  // flatten only
  __m128i indexFlatReg = _mm_set_epi32(intermediates[3][nextIterationToRecalc], 
				       intermediates[2][nextIterationToRecalc], 
				       intermediates[1][nextIterationToRecalc], 
				       intermediates[0][nextIterationToRecalc]);

  // evaluate only
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

  for (size_t i = nextIterationToRecalc; i < dim; i++) {
    //TODO: replace by avx2 gather
    //TODO (better): replace by array access
    __m256d dataTupleReg = _mm256_set_pd(dataTuplePtr[3][i], dataTuplePtr[2][i], 
					 dataTuplePtr[1][i], dataTuplePtr[0][i]);

    //TODO: replace by array access
    __m256d hInverseReg = _mm256_set1_pd((double) hInversePtr[i]);
    __m256d unadjustedReg = _mm256_mul_pd(dataTupleReg, hInverseReg);

    //implies flooring
    __m128i roundedReg = _mm256_cvttpd_epi32(unadjustedReg);	 
    __m128i andedReg = _mm_and_si128(oneIntegerReg, roundedReg);
    __m128i signReg = _mm_xor_si128(oneIntegerReg, andedReg);
    __m128i indexReg = _mm_add_epi32(roundedReg, signReg);	  	  

    //flatten index
    //TODO: can be precomputed as array
    uint32_t actualDirectionGridPoints = hInversePtr[i];
    actualDirectionGridPoints >>= 1;
    __m128i actualDirectionGridPointsReg = _mm_set1_epi32(actualDirectionGridPoints);

    indexFlatReg = _mm_mullo_epi32(indexFlatReg, actualDirectionGridPointsReg);

    __m128i indexShiftedReg = _mm_srli_epi32(indexReg, 1);
	  
    indexFlatReg = _mm_add_epi32(indexFlatReg, indexShiftedReg);
    
    sseUnion.integerRegister = indexFlatReg;
    intermediates[0][i + 1] = sseUnion.uint32Value[0];
    intermediates[1][i + 1] = sseUnion.uint32Value[1];
    intermediates[2][i + 1] = sseUnion.uint32Value[2];
    intermediates[3][i + 1] = sseUnion.uint32Value[3];

    // evaluate
    __m256d indexDoubleReg = _mm256_cvtepi32_pd(indexReg);

    __m256d phi1DEvalReg = _mm256_mul_pd(hInverseReg, dataTupleReg);
    phi1DEvalReg = _mm256_sub_pd(phi1DEvalReg, indexDoubleReg);

    phi1DEvalReg = _mm256_and_pd(phi1DEvalReg, absMask);
    phi1DEvalReg = _mm256_sub_pd(one, phi1DEvalReg);
    phi1DEvalReg = _mm256_max_pd(zero, phi1DEvalReg);

    phiEvalReg = _mm256_mul_pd(phiEvalReg, phi1DEvalReg);
    
    avxUnion.doubleRegister = phiEvalReg;
    evalIndexValues[0][i + 1] = avxUnion.doubleValue[0];
    evalIndexValues[1][i + 1] = avxUnion.doubleValue[1];
    evalIndexValues[2][i + 1] = avxUnion.doubleValue[2];
    evalIndexValues[3][i + 1] = avxUnion.doubleValue[3];

    //may a structure ind[0] im[0] eval[0] null ind[1] ... might help
  }
  _mm_storeu_si128((__m128i *) indexFlat, indexFlatReg);
  _mm256_storeu_pd(phiEval, phiEvalReg);
}
*/

/*
for (size_t i = nextIterationToRecalc; i < dim; i += 4) {	  
  __m256d dataReg = _mm256_loadu_pd(&(dataTuplePtr[i]));
  __m128i hInverseRegInt = _mm_loadu_si128((__m128i *) &(hInversePtr[i]));
  __m256d hInverseReg = _mm256_cvtepi32_pd(hInverseRegInt);
  __m256d unadjustedReg = _mm256_mul_pd(dataReg, hInverseReg);
  __m128i roundedReg = _mm256_cvttpd_epi32(unadjustedReg);	 
  __m128i andedReg = _mm_and_si128(oneIntegerReg, roundedReg);
  __m128i signReg = _mm_xor_si128(oneIntegerReg, andedReg);
  __m128i addedReg = _mm_add_epi32(roundedReg, signReg);
  _mm_storeu_si128((__m128i *) &(indexPtr[i]), addedReg);
 }
*/
