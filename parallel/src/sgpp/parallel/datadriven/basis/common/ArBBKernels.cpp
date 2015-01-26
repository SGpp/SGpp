/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/parallel/datadriven/basis/common/ArBBKernels.hpp>

#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <string>

#include <arbb.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    size_t global_dims;

    arbb::dense<arbb::f32, 2> ArBB_DataSP;
    arbb::dense<arbb::f32, 2> ArBB_LevelSP;
    arbb::dense<arbb::f32, 2> ArBB_IndexSP;

    arbb::dense<arbb::f64, 2> ArBB_Data;
    arbb::dense<arbb::f64, 2> ArBB_Level;
    arbb::dense<arbb::f64, 2> ArBB_Index;

    template <typename fp_Type>
    void arbb_evalGridPoint_oneDim(const arbb::dense<fp_Type>& dataPointsDim, const fp_Type& levelPoint, const fp_Type& index, arbb::dense<fp_Type>& result) {
      result = arbb::max((1.0 - arbb::abs(((levelPoint * dataPointsDim) - index))), 0.0);
    }

    template <typename fp_Type>
    void arbb_evalGridPoint_oneDim_mod(const arbb::dense<fp_Type>& dataPointsDim, const fp_Type& levelPoint, const fp_Type& index, arbb::dense<fp_Type>& result) {
      _if (levelPoint == 2.0) {
        // nothing
      }
      _else_if (index == 1.0) {
        result = arbb::max(2.0 - (levelPoint * dataPointsDim), 0.0);
      }
      _else_if (index == (levelPoint - 1.0)) {
        result = arbb::max((levelPoint * dataPointsDim) - index + 1.0, 0.0);
      }
      _else {
        result = arbb::max((1.0 - arbb::abs(((levelPoint * dataPointsDim) - index))), 0.0);
      } _end_if;
    }

    template <typename fp_Type>
    void arbb_evalTransGridPoint_oneDim(const fp_Type& dataPointDim, const arbb::dense<fp_Type>& level, const arbb::dense<fp_Type>& index, arbb::dense<fp_Type>& result) {
      result = arbb::max((1.0 - arbb::abs(((level * dataPointDim) - index))), 0.0);
    }

    template <typename fp_Type>
    void arbb_evalTransGridPoint_oneDim_mod(const fp_Type& dataPointDim, const arbb::dense<fp_Type>& level, const arbb::dense<fp_Type>& index, arbb::dense<fp_Type>& result) {
      _for(arbb::usize i = 0, i < level.num_rows(), i++) {
        _if (level[i] == 2.0) {
          result[i] = 1.0;
        }
        _else_if (index[i] == 1.0) {
          result[i] = arbb::max(2.0 - (level[i] * dataPointDim), 0.0);
        }
        _else_if (index[i] == (level[i] - 1.0)) {
          result[i] = arbb::max((level[i] * dataPointDim) - index[i] + 1.0, 0.0);
        }
        _else {
          result[i] = arbb::max((1.0 - arbb::abs(((level[i] * dataPointDim) - index[i]))), 0.0);
        } _end_if;
      }
      _end_for;
    }

    template <typename fp_Type>
    void arbb_multTrans(const arbb::dense<fp_Type, 2>& Data, const arbb::dense<fp_Type, 2>& Level, const arbb::dense<fp_Type, 2>& Index, const arbb::dense<fp_Type>& source, arbb::dense<fp_Type>& result) {
      arbb::usize source_size = Data.num_rows();
      arbb::usize storage_size = Level.num_rows();

      arbb::dense<fp_Type, 2> LevelTrans = arbb::transpose(Level);
      arbb::dense<fp_Type, 2> IndexTrans = arbb::transpose(Index);

      _for (arbb::usize i = 0, i < source_size, i++) {
        arbb::dense<fp_Type> curDataPoint = Data.row(i);
        fp_Type s = source[i];

        arbb::dense<fp_Type> temp_result = arbb::fill(s, storage_size);
        arbb::dense<fp_Type> temp;

        // Use normal for loop -> runtime code generation
        for (size_t m = 0; m < global_dims; m++) {
          arbb::usize d = m;
          fp_Type dataDim = curDataPoint[d];

          arbb_evalTransGridPoint_oneDim(dataDim, LevelTrans.row(d), IndexTrans.row(d), temp);

          temp_result *= temp;
        }

        result += temp_result;
      }
      _end_for;
    }

    template <typename fp_Type>
    void arbb_multTrans_mod(const arbb::dense<fp_Type, 2>& Data, const arbb::dense<fp_Type, 2>& Level, const arbb::dense<fp_Type, 2>& Index, const arbb::dense<fp_Type>& source, arbb::dense<fp_Type>& result) {
      arbb::usize source_size = Data.num_rows();
      arbb::usize storage_size = Level.num_rows();

      arbb::dense<fp_Type, 2> LevelTrans = arbb::transpose(Level);
      arbb::dense<fp_Type, 2> IndexTrans = arbb::transpose(Index);

      _for (arbb::usize i = 0, i < source_size, i++) {
        arbb::dense<fp_Type> curDataPoint = Data.row(i);
        fp_Type s = source[i];

        arbb::dense<fp_Type> temp_result = arbb::fill(s, storage_size);
        arbb::dense<fp_Type> temp;

        // Use normal for loop -> runtime code generation
        for (size_t m = 0; m < global_dims; m++) {
          arbb::usize d = m;
          fp_Type dataDim = curDataPoint[d];

          arbb_evalTransGridPoint_oneDim_mod(dataDim, LevelTrans.row(d), IndexTrans.row(d), temp);

          temp_result *= temp;
        }

        result += temp_result;
      }
      _end_for;
    }

    template <typename fp_Type>
    void arbb_mult(const arbb::dense<fp_Type, 2>& Data, const arbb::dense<fp_Type, 2>& Level, const arbb::dense<fp_Type, 2>& Index, const arbb::dense<fp_Type>& alpha, arbb::dense<fp_Type>& result) {
      arbb::usize result_size = result.length();
      arbb::usize storage_size = Level.num_rows();

      arbb::dense<fp_Type, 2> DataTrans = arbb::transpose(Data);

      _for (arbb::usize j = 0, j < storage_size, j++) {
        arbb::dense<fp_Type> curLevel = Level.row(j);
        arbb::dense<fp_Type> curIndex = Index.row(j);
        fp_Type a = alpha[j];

        arbb::dense<fp_Type> temp_result = arbb::fill(a, result_size);
        arbb::dense<fp_Type> temp;

        // Use normal for loop -> runtime code generation
        for (size_t m = 0; m < global_dims; m++) {
          arbb::usize d = m;
          fp_Type l = curLevel[d];
          fp_Type i = curIndex[d];

          //arbb::dense<fp_Type> index = arbb::fill(i, result_size);

          arbb_evalGridPoint_oneDim(DataTrans.row(d), l, i, temp);

          temp_result *= temp;
        }

        result += temp_result;
      }
      _end_for;
    }

    template <typename fp_Type>
    void arbb_mult_mod(const arbb::dense<fp_Type, 2>& Data, const arbb::dense<fp_Type, 2>& Level, const arbb::dense<fp_Type, 2>& Index, const arbb::dense<fp_Type>& alpha, arbb::dense<fp_Type>& result) {
      arbb::usize result_size = result.length();
      arbb::usize storage_size = Level.num_rows();

      arbb::dense<fp_Type, 2> DataTrans = arbb::transpose(Data);

      _for (arbb::usize j = 0, j < storage_size, j++) {
        arbb::dense<fp_Type> curLevel = Level.row(j);
        arbb::dense<fp_Type> curIndex = Index.row(j);
        fp_Type a = alpha[j];

        arbb::dense<fp_Type> temp_result = arbb::fill(a, result_size);
        arbb::dense<fp_Type> temp = arbb::fill(1.0, result_size);

        // Use normal for loop -> runtime code generation
        for (size_t m = 0; m < global_dims; m++) {
          arbb::usize d = m;
          fp_Type l = curLevel[d];
          fp_Type i = curIndex[d];

          //arbb::dense<fp_Type> index = arbb::fill(i, result_size);

          arbb_evalGridPoint_oneDim_mod(DataTrans.row(d), l, i, temp);

          temp_result *= temp;
        }

        result += temp_result;
      }
      _end_for;
    }

    ArBBKernels::ArBBKernels() {
      isMultTransSPfirst = true;
      isMultSPfirst = true;

      isMultTransfirst = true;
      isMultfirst = true;
    }

    ArBBKernels::~ArBBKernels() {
    }

    double ArBBKernels::multTransArBB(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims = dims;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_source;

        if (isMultTransfirst && isMultfirst) {
          arbb::bind(ArBB_Data, ptrData, dims, sourceSize);
          arbb::bind(ArBB_Level, ptrLevel, dims, storageSize);
          arbb::bind(ArBB_Index, ptrIndex, dims, storageSize);
          isMultTransfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans<arbb::f64>))(ArBB_Data, ArBB_Level, ArBB_Index, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels::multArBB(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims = dims;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_alpha;

        if (isMultTransfirst && isMultfirst) {
          arbb::bind(ArBB_Data, ptrData, dims, result_size);
          arbb::bind(ArBB_Level, ptrLevel, dims, storageSize);
          arbb::bind(ArBB_Index, ptrIndex, dims, storageSize);
          isMultfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult<arbb::f64>))(ArBB_Data, ArBB_Level, ArBB_Index, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels::multTransSPArBB(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims = dims;

      try {
        arbb::dense<arbb::f32> ArBB_result;
        arbb::dense<arbb::f32> ArBB_source;

        if (isMultTransSPfirst && isMultSPfirst) {
          arbb::bind(ArBB_DataSP, ptrData, dims, sourceSize);
          arbb::bind(ArBB_LevelSP, ptrLevel, dims, storageSize);
          arbb::bind(ArBB_IndexSP, ptrIndex, dims, storageSize);
          isMultTransSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans<arbb::f32>))(ArBB_DataSP, ArBB_LevelSP, ArBB_IndexSP, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels::multSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims = dims;

      try {
        arbb::dense<arbb::f32> ArBB_result;
        arbb::dense<arbb::f32> ArBB_alpha;

        if (isMultTransSPfirst && isMultSPfirst) {
          arbb::bind(ArBB_DataSP, ptrData, dims, result_size);
          arbb::bind(ArBB_LevelSP, ptrLevel, dims, storageSize);
          arbb::bind(ArBB_IndexSP, ptrIndex, dims, storageSize);
          isMultSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult<arbb::f32>))(ArBB_DataSP, ArBB_LevelSP, ArBB_IndexSP, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels::multModTransArBB(double* ptrSource, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims = dims;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_source;

        if (isMultModTransfirst && isMultModfirst) {
          arbb::bind(ArBB_Data, ptrData, dims, sourceSize);
          arbb::bind(ArBB_Level, ptrLevel, dims, storageSize);
          arbb::bind(ArBB_Index, ptrIndex, dims, storageSize);
          isMultModTransfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans_mod<arbb::f64>))(ArBB_Data, ArBB_Level, ArBB_Index, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels::multModArBB(double* ptrAlpha, double* ptrData, double* ptrLevel, double* ptrIndex, double* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims = dims;

      try {
        arbb::dense<arbb::f64> ArBB_result;
        arbb::dense<arbb::f64> ArBB_alpha;

        if (isMultModTransfirst && isMultModfirst) {
          arbb::bind(ArBB_Data, ptrData, dims, result_size);
          arbb::bind(ArBB_Level, ptrLevel, dims, storageSize);
          arbb::bind(ArBB_Index, ptrIndex, dims, storageSize);
          isMultModfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult_mod<arbb::f64>))(ArBB_Data, ArBB_Level, ArBB_Index, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels::multModTransSPArBB(float* ptrSource, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrGlobalResult, size_t sourceSize, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims = dims;

      try {
        arbb::dense<arbb::f32> ArBB_result;
        arbb::dense<arbb::f32> ArBB_source;

        if (isMultModTransSPfirst && isMultModSPfirst) {
          arbb::bind(ArBB_DataSP, ptrData, dims, sourceSize);
          arbb::bind(ArBB_LevelSP, ptrLevel, dims, storageSize);
          arbb::bind(ArBB_IndexSP, ptrIndex, dims, storageSize);
          isMultModTransSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrGlobalResult, storageSize);
        arbb::bind(ArBB_source, ptrSource, sourceSize);

        arbb::call(&(arbb_multTrans_mod<arbb::f32>))(ArBB_DataSP, ArBB_LevelSP, ArBB_IndexSP, ArBB_source, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    double ArBBKernels::multModSPArBB(float* ptrAlpha, float* ptrData, float* ptrLevel, float* ptrIndex, float* ptrResult, size_t result_size, size_t storageSize, size_t dims) {
      double time = 0.0;

      global_dims = dims;

      try {
        arbb::dense<arbb::f32> ArBB_result;
        arbb::dense<arbb::f32> ArBB_alpha;

        if (isMultModTransSPfirst && isMultModSPfirst) {
          arbb::bind(ArBB_DataSP, ptrData, dims, result_size);
          arbb::bind(ArBB_LevelSP, ptrLevel, dims, storageSize);
          arbb::bind(ArBB_IndexSP, ptrIndex, dims, storageSize);
          isMultModSPfirst = false;
        }

        arbb::bind(ArBB_result, ptrResult, result_size);
        arbb::bind(ArBB_alpha, ptrAlpha, storageSize);

        arbb::call(&(arbb_mult_mod<arbb::f32>))(ArBB_DataSP, ArBB_LevelSP, ArBB_IndexSP, ArBB_alpha, ArBB_result);
      } catch (const std::exception& e) {
        std::cout << "Error using Intel ArBB: " << e.what() << std::endl;
      }

      return time;
    }

    void ArBBKernels::resetKernels() {
      isMultTransSPfirst = true;
      isMultSPfirst = true;
      isMultModTransSPfirst = true;
      isMultModSPfirst = true;

      isMultTransfirst = true;
      isMultfirst = true;
      isMultModTransfirst = true;
      isMultModfirst = true;
    }

  }

}
