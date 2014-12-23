static inline void flattenIndexVectorized(uint32_t *(&intermediates)[4], size_t dim, 
					  uint32_t *maxIndicesPtr, uint32_t *(&indexPtr)[4], 
					  size_t toRecalc, uint32_t (&indexFlat)[4]) {

  union {
    __m128d doubleRegister;
    __m128i integerRegister;
    uint32_t uint32Value[4];
  } sseUnion;
	
  __m128i indexFlatReg = _mm_set_epi32(intermediates[3][toRecalc], intermediates[2][toRecalc], 
				       intermediates[1][toRecalc], intermediates[0][toRecalc]);

  //cout << "<<<<<<<<<<ref<<<<<<<<<<" << endl;
	
  for (size_t i = toRecalc; i < dim; i++) {
    uint32_t actualDirectionGridPoints = maxIndicesPtr[i];
    actualDirectionGridPoints >>= 1;

    __m128i actualDirectionGridPointsReg = _mm_set1_epi32(actualDirectionGridPoints);

    // sseUnion.integerRegister = indexFlatReg;
    // cout << "combFlat[3] = " << sseUnion.uint32Value[3] << endl;

    indexFlatReg = _mm_mullo_epi32(indexFlatReg, actualDirectionGridPointsReg);

    // sseUnion.integerRegister = indexFlatReg;
    // cout << "multPost[3] = " << sseUnion.uint32Value[3] << endl;

    __m128i actualIndexReg = _mm_set_epi32(indexPtr[3][i], indexPtr[2][i], 
					   indexPtr[1][i], indexPtr[0][i]);

    // sseUnion.integerRegister = actualIndexReg;
    // cout << "combFlat[3] = " << sseUnion.uint32Value[3] << endl;

    actualIndexReg = _mm_srli_epi32(actualIndexReg, 1);
	  
    indexFlatReg = _mm_add_epi32(indexFlatReg, actualIndexReg);
    
    sseUnion.integerRegister = indexFlatReg;
    intermediates[0][i + 1] = sseUnion.uint32Value[0];
    intermediates[1][i + 1] = sseUnion.uint32Value[1];
    intermediates[2][i + 1] = sseUnion.uint32Value[2];
    intermediates[3][i + 1] = sseUnion.uint32Value[3];
    //cout << "combFlat[3] = " << sseUnion.uint32Value[3] << endl;
    
  }

  _mm_storeu_si128((__m128i *) indexFlat, indexFlatReg);
}
